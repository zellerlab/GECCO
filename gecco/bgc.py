import csv
import typing
from typing import List, Optional

import numpy as np
from gecco.preprocessing import flatten


class Protein(object):
    """A single protein withing a BGC.
    """

    def __init__(
        self,
        seq_id: str,
        start: int,
        end: int,
        name: str,
        probability: float = 0.0,
        domains: Optional[List[str]] = None,
        weights: Optional[List[float]] = None,
    ) -> None:
        """Create a new `Protein` object.

        For practical reasons the strandedness of the protein is not taken
        into account, i.e. the ``start`` attribute will be the lowest of the
        two indices delimiting the protein.

        Arguments:
            seq_id (str): The identifier of the sequence containing the
                protein.
            start (int): The location where the protein starts.
            end (int): The location where the protein ends.
            name (str): The identifier of the protein itself.
            probability (float): The probability with which the protein was
                predicted.
            domains (iterable of `str`): The domains of the protein, if any
                were predicted.
            weights (iterable of `float`): The prediction weights associated
                with each domain of the protein.

        Raises:
            ValueError: when `domains` and `weights` are not the same size.
        """
        self.seq_id = seq_id
        self.start = min(start, end)
        self.end = max(start, end)
        self.name = name
        self.probability: float = probability

        self.domains = np.array([] if domains is None else list(domains))
        if weights is not None:
            self.weights: np.ndarray = np.array(list(weights))
        else:
            self.weights = np.ones(len(self.domains))
        if len(self.weights) != len(self.domains):
            raise ValueError("length of `domains` and `weights` differs")

    def is_potential_cluster(self, thresh: float = 0.3) -> bool:
        """Check whether or not this protein is a potential cluster.

        Arguments:
            thresh (float): The probability threshold above which a protein is
                considered a potential cluster. Defaults to `0.3`.
        """
        return self.probability > thresh


class BGC(object):
    """A biosynthetic gene cluster, containing multiple proteins.
    """

    def __init__(
        self,
        proteins: List[Protein],
        name: Optional[str] = None,
        bgc_type: Optional[str] = None,
        type_prob: Optional[float] = None,
    ) -> None:
        """Create a `BGC` from a list of `Protein` objects.

        Arguments:
            proteins (`list` of `Protein`): The proteins contained in the
                current BGC.
            name (`str`, optional): A name to use to identify the BGC. If none
                is given, then the sequence id of any `Protein` will be used.
            bgc_type (`str`, optional): The BGC type, if it is known.
            type_prob (`float`, optional): The probability with which the BGC
                type was predicted, if any.

        """
        self.proteins = proteins
        self.seq_id = self.proteins[0].seq_id
        self.prot_ids = np.array([p.name for p in proteins])
        self.start = min(p.start for p in proteins)
        self.end = max(p.end for p in proteins)
        self.domains = np.array([p.domains for p in proteins])
        self.weights = np.array([p.weights for p in proteins])
        self.probabilities = np.array([p.probability for p in proteins])
        self.name = name if name else self.seq_id
        self.bgc_type = bgc_type
        self.type_prob = type_prob

    def is_valid(self, criterion: str = "antismash", thresh: float = 0.6) -> bool:
        if criterion == "antismash":
            # These are the default options only
            return self._antismash_check(n_biopfams=5, p_thresh=0.6, n_cds=5)
        elif criterion == "gecco":
            return self._gecco_check()
        else:
            raise ValueError(f"invalid criterion: {criterion!r}")

    def write_to_file(self, handle: typing.TextIO, long: bool = False) -> None:
        """Write the BGC as a single row to a text file.

        Arguments:
            handle (file handle): A file handle open in text mode.
            long (bool, optional): Whether or not to include additional
                information (such as the ids of the proteins or of the domains
                within the BGC). Defaults to `False`.
        """
        probs = np.hstack(self.probabilities)
        row = [
            self.seq_id,
            self.name,
            self.start,
            self.end,
            probs.mean(),
            probs.max(),
            self.bgc_type,
            self.type_prob,
        ]
        if long:
            row.append(",".join(np.hstack(self.prot_ids)))  # prots
            row.append(",".join(np.hstack(self.domains)))  # pfam
        row_str = [str(item) for item in row]
        csv.writer(handle, dialect="excel-tab").writerow(row_str)

    def domain_composition(self, all_possible: Optional[np.ndarray] = None) -> np.ndarray:
        """Computes weighted domain composition with respect to ``all_possible``.

        Arguments:
            all_possible (`~numpy.ndarray`): An array containing all domain
                names to consider when computing domain composition for the
                BGC. If `None` given, then all domains are taken into account.

        Returns:
            `~numpy.ndarray`: A numerical array containing the relative domain
            composition of the BGC.
        """
        doms = np.hstack(self.domains)
        w = np.hstack(self.weights)
        if all_possible is None:
            all_possible = np.unique(doms)
        comp_arr = np.zeros(len(all_possible))
        for i, dom in enumerate(all_possible):
            n = list(doms).count(dom)
            weight = w[doms == dom].mean() if n > 0 else 0
            comp_arr[i] = n * weight
        comp_arr = comp_arr / comp_arr.sum()
        return np.nan_to_num(comp_arr, copy=False)

    def domain_counts(self, all_possible=None) -> np.ndarray:
        """Computes domain counts with respect to ``all_possible``.
        """
        doms = list(np.hstack(self.domains))
        if all_possible is None:
            all_possible = np.unique(doms)
        return np.array([doms.count(d) for d in all_possible])

    def _antismash_check(self, n_biopfams=5, p_thresh=0.6, n_cds=5) -> bool:
        """
        Checks for cluster validity using antiSMASH criteria:
        1) MEAN p over thresh
        2) At least n_biopfams intersection between bio_pfams and the pfams in the
        cluster
        """

        # fmt: off
        bio_pfams = {
            "PF00109", "PF02801", "PF08659", "PF00378", "PF08541", "PF08545",
            "PF02803", "PF00108", "PF02706", "PF03364", "PF08990", "PF00501",
            "PF00668", "PF08415", "PF00975", "PF03061", "PF00432", "PF00494",
            "PF03936", "PF01397", "PF00432", "PF04275", "PF00348", "PF02401",
            "PF04551", "PF00368", "PF00534", "PF00535", "PF02922", "PF01041",
            "PF00128", "PF00908", "PF02719", "PF04321", "PF01943", "PF02806",
            "PF02350", "PF02397", "PF04932", "PF01075", "PF00953", "PF01050",
            "PF03033", "PF01501", "PF05159", "PF04101", "PF02563", "PF08437",
            "PF02585", "PF01721", "PF02052", "PF02674", "PF03515", "PF04369",
            "PF08109", "PF08129", "PF09221", "PF09683", "PF10439", "PF11420",
            "PF11632", "PF11758", "PF12173", "PF04738", "PF04737", "PF04604",
            "PF05147", "PF08109", "PF08129", "PF08130", "PF00155", "PF00202",
            "PF00702", "PF06339", "PF04183", "PF10331", "PF03756", "PF00106",
            "PF01370", "PF00107", "PF08240", "PF00441", "PF02770", "PF02771",
            "PF08028", "PF01408", "PF02894", "PF00984", "PF00725", "PF03720",
            "PF03721", "PF07993", "PF02737", "PF00903", "PF00037", "PF04055",
            "PF00171", "PF00067", "PF01266", "PF01118", "PF02668", "PF00248",
            "PF01494", "PF01593", "PF03992", "PF00355", "PF01243", "PF00384",
            "PF01488", "PF00857", "PF04879", "PF08241", "PF08242", "PF00698",
            "PF00483", "PF00561", "PF00583", "PF01636", "PF01039", "PF00288",
            "PF00289", "PF02786", "PF01757", "PF02785", "PF02409", "PF01553",
            "PF02348", "PF00891", "PF01596", "PF04820", "PF02522", "PF08484",
            "PF08421",
        }

        bio_crit = len(set(np.hstack(self.domains)) & bio_pfams) >= n_biopfams
        p_crit = np.mean([p.mean() for p in self.probabilities]) >= p_thresh
        cds_crit = len(np.hstack(self.prot_ids)) >= n_cds

        return bio_crit and p_crit and cds_crit

    def _gecco_check(self) -> bool:
        return True
