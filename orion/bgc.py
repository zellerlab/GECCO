import numpy as np
from orion.preprocessing import flatten

class Protein(object):
    """Definition of a Protein"""
    def __init__(self, start, end, name, p=0.0, domains=[], weights=None):
        self.start = min(start, end)
        self.end = max(start, end)
        self.name = np.array(name)
        self.domains = np.array(domains)
        if weights is not None:
            self.weights = np.array(weights)
        else:
            self.weights = np.ones(len(domains))
        self.probs = np.array(p)

    def is_potential_cluster(self, thresh=0.3):
        return self.probs.mean() > thresh


class BGC(object):
    """A biosynthetic gene cluster with multiple proteins"""
    def __init__(self, proteins, name="cluster", bgc_type=None):
        self.name = name
        self.type = bgc_type
        self.proteins = proteins
        self.prot_ids = np.array([p.name for p in proteins])
        self.start = min([p.start for p in proteins])
        self.end = max([p.end for p in proteins])
        self.domains = np.array([p.domains for p in proteins])
        self.weights = np.array([p.weights for p in proteins])
        self.probs = np.array([p.probs for p in proteins])

    def is_valid(self, criterion="antismash", thresh=0.6):
        if criterion == "antismash":
            # These are the default options only
            return self._antismash_check(
                n_biopfams = 5,
                p_thresh = 0.6,
                n_cds = 5
            )
        if criterion == "orion":
            return self._orion_check()

    def write_to_file(self, handle, short=True):
        p_mean = np.hstack(self.probs).mean()
        p_max = np.hstack(self.probs).max()
        if short:
            row = [self.name, self.start, self.end, p_mean, p_max, self.type]
        else:
            prot = ",".join(np.hstack(self.prot_ids))
            pfam = ",".join(np.hstack(self.domains))
            row = [self.name, self.start, self.end, p_mean, p_max, self.type, prot, pfam]
        row = map(str, row)
        handle.write("\t".join(row) + "\n")

    def domain_composition(self, all_possible=None):
        """Computes weighted domain composition with respect to all_possible.
        """
        doms = np.hstack(self.domains)
        w = np.hstack(self.weights)
        if all_possible is None:
            all_possible = np.unique(doms)
        comp_arr = np.zeros(len(all_possible))
        for i in range(len(all_possible)):
            n = list(doms).count(all_possible[i])
            if n > 0:
                weight = w[doms == all_possible[i]].mean()
            else:
                weight = 0
            comp_arr[i] = n * weight
        return comp_arr / comp_arr.sum()

    def domain_counts(self, all_possible=None):
        """Comoutes domain counts with respect to all_possible.
        """
        doms = np.hstack(self.domains)
        if all_possible is None:
            all_possible = np.unique(doms)
        comp_arr = np.zeros(len(all_possible))
        for i in range(len(all_possible)):
            n = list(doms).count(all_possible[i])
            comp_arr[i] = n
        return comp_arr

    def _antismash_check(self, n_biopfams=5, p_thresh=0.6, n_cds=5):
        """
        Checks for cluster validity using antiSMASH criteria:
        1) MEAN p over thresh
        2) At least n_biopfams intersection between bio_pfams and the pfams in the
        cluster
        """

        bio_pfams = {
            "PF00109", "PF02801", "PF08659", "PF00378", "PF08541",
            "PF08545", "PF02803", "PF00108", "PF02706", "PF03364", "PF08990", "PF00501",
            "PF00668", "PF08415", "PF00975", "PF03061", "PF00432", "PF00494", "PF03936",
            "PF01397", "PF00432", "PF04275", "PF00348", "PF02401", "PF04551", "PF00368",
            "PF00534", "PF00535", "PF02922", "PF01041","PF00128", "PF00908","PF02719", "PF04321", "PF01943", "PF02806", "PF02350", "PF02397", "PF04932","PF01075",
            "PF00953","PF01050", "PF03033", "PF01501", "PF05159", "PF04101", "PF02563",
            "PF08437", "PF02585", "PF01721", "PF02052", "PF02674","PF03515", "PF04369",
            "PF08109", "PF08129", "PF09221", "PF09683", "PF10439", "PF11420", "PF11632",
            "PF11758", "PF12173","PF04738", "PF04737", "PF04604", "PF05147", "PF08109",
            "PF08129", "PF08130", "PF00155", "PF00202", "PF00702", "PF06339","PF04183",
            "PF10331", "PF03756", "PF00106", "PF01370", "PF00107", "PF08240", "PF00441",
            "PF02770", "PF02771", "PF08028","PF01408", "PF02894", "PF00984", "PF00725",
            "PF03720", "PF03721", "PF07993", "PF02737", "PF00903", "PF00037", "PF04055",
            "PF00171", "PF00067", "PF01266", "PF01118", "PF02668", "PF00248", "PF01494",
            "PF01593", "PF03992", "PF00355", "PF01243","PF00384", "PF01488", "PF00857",
            "PF04879", "PF08241", "PF08242", "PF00698", "PF00483", "PF00561", "PF00583",
            "PF01636","PF01039", "PF00288", "PF00289", "PF02786", "PF01757", "PF02785",
            "PF02409", "PF01553", "PF02348", "PF00891", "PF01596","PF04820", "PF02522",
            "PF08484", "PF08421"
        }

        bio_crit = len(set(np.hstack(self.domains)) & bio_pfams) >= n_biopfams
        p_crit = np.mean([p.mean() for p in self.probs]) >= p_thresh
        cds_crit = len(np.hstack(self.prot_ids)) >= n_cds

        return (bio_crit and p_crit and cds_crit)

    def _orion_check(self):
        return True
