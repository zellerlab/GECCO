"""Implementation of the ``gecco train`` subcommand.
"""

import csv
import hashlib
import itertools
import logging
import multiprocessing.pool
import os
import operator
import pickle
import random
import typing

import numpy
import scipy.sparse
import tqdm
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from ._base import Command
from ...model import FeatureTable, ClusterTable, Cluster
from ...crf import ClusterCRF
from ...refine import ClusterRefiner
from ...hmmer import HMMER


class Train(Command):  # noqa: D101

    summary = "train the CRF model on an embedded feature table."

    @classmethod
    def doc(cls, fast=False):
        return f"""
        gecco train - {cls.summary}

        Usage:
            gecco train (-h | --help)
            gecco train --features <table> --clusters <table> [options]

        Arguments:
            -f <data>, --features <table>   a domain annotation table, used to
                                            train the CRF model.
            -c <data>, --clusters <table>   a cluster annotation table, used to
                                            extract the domain composition for
                                            the type classifier.

        Parameters:
            -o <out>, --output-dir <out>    the directory to use for the model
                                            files. [default: CRF]
            -j <jobs>, --jobs <jobs>        the number of CPUs to use for
                                            multithreading. Use 0 to use all of the
                                            available CPUs. [default: 0]

        Parameters - Domain Annotation:
            -e <e>, --e-filter <e>          the e-value cutoff for domains to
                                            be included [default: 1e-5]

        Parameters - Training:
            --c1 <C1>                       parameter for L1 regularisation.
                                            [default: 0.15]
            --c2 <C2>                       parameter for L2 regularisation.
                                            [default: 0.15]
            --feature-type <type>           how features should be extracted
                                            (single, overlap, or group).
                                            [default: group]
            --truncate <N>                  the maximum number of rows to use from
                                            the training set.
            --overlap <N>                   how much overlap to consider if
                                            features overlap. [default: 2]
            --no-shuffle                    disable shuffling of the data before
                                            fitting the model.
            --select <N>                    fraction of most significant features
                                            to select from the training data.

        """

    def _check(self) -> typing.Optional[int]:
        retcode = super()._check()
        if retcode is not None:
            return retcode

        # Check the input exists
        for input_ in self.args["--features"], self.args["--clusters"]:
            if not os.path.exists(input_):
                self.logger.error("could not locate input file: {!r}", input_)
                return 1

        # Check the `--feature-type`
        type_ = self.args["--feature-type"]
        if type_ not in {"single", "overlap", "group"}:
            self.logger.error("Invalid value for `--feature-type`: {}", type_)
            return 1

        # Check value of numeric arguments
        if self.args["--truncate"] is not None:
            self.args["--truncate"] = int(self.args["--truncate"])
        self.args["--overlap"] = int(self.args["--overlap"])
        self.args["--c1"] = float(self.args["--c1"])
        self.args["--c2"] = float(self.args["--c2"])
        self.args["--e-filter"] = e_filter = float(self.args["--e-filter"])
        if e_filter < 0 or e_filter > 1:
            self.logger.error("Invalid value for `--e-filter`: {}", e_filter)
            return 1
        if self.args["--select"] is not None:
            self.args["--select"] = float(self.args["--select"])

        # Check the `--jobs`flag
        self.args["--jobs"] = jobs = int(self.args["--jobs"])
        if jobs == 0:
            self.args["--jobs"] = multiprocessing.cpu_count()

        return None

    def __call__(self) -> int:  # noqa: D102
        # --- LOADING AND PREPROCESSING --------------------------------------
        # Load the table
        self.logger.info("Loading the data")
        with open(self.args["--features"]) as in_:
            features = FeatureTable.load(in_)

        # Converting table to genes and sort by location
        self.logger.info("Sorting genes by location")
        genes = sorted(features.to_genes(), key=operator.attrgetter("source.id", "start", "end"))
        for gene in genes:
            gene.protein.domains.sort(key=operator.attrgetter("start", "end"))
        del features

        # --- MODEL FITTING --------------------------------------------------
        self.logger.info("Fitting the CRF model to the training data")
        crf = ClusterCRF(
            self.args["--feature-type"],
            algorithm="lbfgs",
            overlap=self.args["--overlap"],
            c1=self.args["--c1"],
            c2=self.args["--c2"],
        )
        crf.fit(genes, select=self.args["--select"], shuffle=not self.args["--no-shuffle"])

        # --- MODEL SAVING ---------------------------------------------------
        os.makedirs(self.args["--output-dir"], exist_ok=True)
        model_out = os.path.join(self.args["--output-dir"], "model.pkl")
        self.logger.info("Writing the model to {!r}", model_out)
        with open(model_out, "wb") as out:
            pickle.dump(crf, out, protocol=4)

        self.logger.debug("Computing and saving model checksum")
        hasher = hashlib.md5()
        with open(model_out, "rb") as out:
            hasher.update(out.read())  # FIXME: iterate on file blocks
        with open(f"{model_out}.md5", "w") as out_hash:
            out_hash.write(hasher.hexdigest())

        self.logger.info("Writing transitions weights")
        with open(os.path.join(self.args["--output-dir"], "model.trans.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["from", "to", "weight"])
            for labels, weight in crf.model.transition_features_.items():
                writer.writerow([*labels, weight])

        self.logger.info("Writing state weights")
        with open(os.path.join(self.args["--output-dir"], "model.state.tsv"), "w") as f:
            writer = csv.writer(f, dialect="excel-tab")
            writer.writerow(["attr", "label", "weight"])
            for attrs, weight in crf.model.state_features_.items():
                writer.writerow([*attrs, weight])

        # --- DOMAIN COMPOSITION ---------------------------------------------
        self.logger.info("Extracting domain composition for type classifier")

        self.logger.debug("Loading cluster table")
        with open(self.args["--clusters"]) as f:
            clusters = [
                Cluster(c.bgc_id, [g for g in genes if g.id in c.proteins], c.type)
                for c in sorted(ClusterTable.load(f), key=operator.attrgetter("bgc_id"))
            ]

        self.logger.debug("Finding the array of possible domains")
        if crf.significant_features is not None:
            all_possible = sorted(crf.significant_features)
        else:
            all_possible = sorted({d.name for g in genes for d in g.protein.domains})

        self.logger.debug("Saving training matrix labels for BGC type classifier")
        with open(os.path.join(self.args["--output-dir"], "domains.tsv"), "w") as out:
            out.writelines([f"{domain}\n" for domain in all_possible])
        with open(os.path.join(self.args["--output-dir"], "types.tsv"), "w") as out:
            writer = csv.writer(out, dialect="excel-tab")
            for cluster in clusters:
                types =  ";".join(ty.name for ty in cluster.type.unpack())
                writer.writerow([cluster.id, types])

        self.logger.debug("Building new domain composition matrix")
        comp = numpy.array([c.domain_composition(all_possible) for c in clusters])

        comp_out = os.path.join(self.args["--output-dir"], "compositions.npz")
        self.logger.debug("Saving new domain composition matrix to {!r}", comp_out)
        scipy.sparse.save_npz(comp_out, scipy.sparse.coo_matrix(comp))

        return 0
