#!/usr/bin/env python

import configparser
import csv
import glob
import gzip
import hashlib
import io
import json
import math
import os
import re
import shutil
import ssl
import sys
import tarfile
import urllib.request
from functools import partial
from xml.etree import ElementTree as etree

import setuptools
from distutils.command.build import build as _build
from setuptools.command.sdist import sdist as _sdist
from tqdm import tqdm
from pyhmmer.plan7 import HMMFile


class sdist(_sdist):
    """An extension to the `sdist` command that generates a `pyproject.toml`.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", "build-backend", '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)

        # run the rest of the packaging
        _sdist.run(self)


class update_model(setuptools.Command):
    """A custom command to update the internal CRF model.
    """

    description = "update the CRF model embedded in the source"
    user_options = [
        ("model=", "m", "the path to the new CRF model to use"),
    ]

    def initialize_options(self):
        self.model = None
        self.domain = None

    def finalize_options(self):
        if self.model is None:
            raise ValueError("--model argument must be given")
        elif not os.path.isdir(self.model):
            raise FileNotFoundError(self.model)

    def info(self, msg):
        self.announce(msg, level=2)

    def run(self):
        # Copy the file to the new in-source location and compute its hash.
        hasher = hashlib.md5()
        self.info("Copying the trained CRF model to the in-source location")
        with open(os.path.join(self.model, "model.pkl"), "rb") as src:
            with open(os.path.join("gecco", "crf", "model.pkl"), "wb") as dst:
                read = lambda: src.read(io.DEFAULT_BUFFER_SIZE)
                for chunk in iter(read, b""):
                    hasher.update(chunk)
                    dst.write(chunk)

        # Write the hash to the signature file next to the model
        self.info("Writing the MD5 signature file")
        with open(os.path.join("gecco", "crf", "model.pkl.md5"), "w") as sig:
            sig.write(hasher.hexdigest())

        # Update the domain composition table
        self.info("Copying the KNN training data to the in-source location")
        for filename in ["compositions.npz", "domains.tsv", "types.tsv"]:
            src = os.path.join(self.model, filename)
            dst = os.path.join("gecco", "knn", filename)
            shutil.copy(src=src, dst=dst)

        # Update the interpro entries
        path = os.path.join("gecco", "interpro", "interpro.json.gz")
        self.info("getting Pfam entries from InterPro")
        entries = self.download_interpro_entries("pfam")
        self.info("getting Tigrfam entries from InterPro")
        entries.extend(self.download_interpro_entries("tigrfams"))
        with gzip.open(path, "wt") as dest:
            json.dump(entries, dest)

    def download_interpro_entries(self, db):
        next = "https://www.ebi.ac.uk:443/interpro/api/entry/all/{}/?page_size=200".format(db)
        entries = []
        context = ssl._create_unverified_context()
        with tqdm(desc=db, leave=False) as pbar:
            while next:
                with urllib.request.urlopen(next, context=context) as res:
                    payload = json.load(res)
                    pbar.total = payload["count"]
                    next = payload["next"]
                    entries.extend(payload["results"])
                    pbar.update(len(payload["results"]))
        return entries


class build_data(setuptools.Command):
    """A custom `setuptools` command to download data before wheel creation.
    """

    description = "download the HMM libraries used by GECCO to annotate proteins"
    user_options = [
        ("inplace", "i", "ignore build-lib and put data alongside your Python code")
    ]

    def initialize_options(self):
        self.inplace = False

    def finalize_options(self):
        _build_py = self.get_finalized_command("build_py")
        self.build_lib = _build_py.build_lib

    def info(self, msg):
        self.announce(msg, level=2)

    def run(self):
        self.mkpath(self.build_lib)

        domains_file = os.path.join("gecco", "types", "domains.tsv")
        self.info("loading domain accesssions from {}".format(domains_file))
        with open(domains_file, "rb") as f:
            domains = [line.strip() for line in f]

        for in_ in glob.iglob(os.path.join("gecco", "hmmer", "*.ini")):
            local = os.path.join(self.build_lib, in_).replace(".ini", ".h3m")
            self.mkpath(os.path.dirname(local))
            self.make_file([in_], local, self.download, (in_, domains))
            if self.inplace:
                copy = in_.replace(".ini", ".h3m")
                self.make_file([local], copy, shutil.copy, (local, copy))

    def download(self, in_, domains):
        cfg = configparser.ConfigParser()
        cfg.read(in_)
        out = os.path.join(self.build_lib, in_.replace(".ini", ".h3m"))

        try:
            self.download_hmm(out, domains, dict(cfg.items("hmm")))
        except:
            if os.path.exists(out):
                os.remove(out)
            raise

    def download_hmm(self, output, domains, options):
        base = "https://github.com/zellerlab/GECCO/releases/download/v{version}/{id}.hmm.gz"
        url = base.format(id=options["id"], version=self.distribution.get_version())
        # attempt to use the GitHub releases URL, otherwise fallback to official URL
        try:
            self.announce("fetching {}".format(url), level=2)
            response = urllib.request.urlopen(url)
        except urllib.error.HTTPError:
            self.announce("using fallback {}".format(options["url"]), level=2)
            response = urllib.request.urlopen(options["url"])
        # download the HMM
        format = dict(
            total=int(response.headers["Content-Length"]),
            desc=os.path.basename(output),
            leave=False,
        )
        with tqdm.wrapattr(response, "read", **format) as src:
            with open(output, "wb") as dst:
                nwritten = 0
                for hmm in HMMFile(gzip.open(src)):
                    if hmm.accession.split(b".")[0] in domains:
                        hmm.write(dst, binary=True)
                        nwritten += 1


class build(_build):
    """A hacked `build` command that will also run `build_data`.
    """

    def run(self):
        self.run_command("build_data")
        _build.run(self)


if __name__ == "__main__":
    setuptools.setup(
        cmdclass={
            "build": build,
            "build_data": build_data,
            "sdist": sdist,
            "update_model": update_model,
        },
    )
