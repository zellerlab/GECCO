import pandas as pd
import numpy as np

class ClusterRefiner(object):

    def __init__(self, threshold=0.5, biosynthetic_pfams=5, seq_col="sequence_id",
            prot_col="protein_id", min_domains=1, min_proteins=5, join_width=1):

        self.thresh = threshold
        self.n_biopfams = biosynthetic_pfams
        self.n_domains = min_domains
        self.n_proteins = min_proteins
        self.n_proteins = min_proteins
        self.join_width = join_width
        self.seq_col = seq_col
        self.prot_col = prot_col
        self.grouping = [seq_col, prot_col]

    def refine(self, pfam_df, target_col="Y_pred", method="antismash"):
        if method == "antismash":
            clusters = self._antismash_refine(pfam_df, target_col)
            return clusters

    def _antismash_refine(self, pfam_df, target_col):

        bio_pfams = {"PF00109", "PF02801", "PF08659", "PF00378", "PF08541",
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
        "PF08484", "PF08421"}

        cluster_condition = {
            target_col: lambda r: 1 if (any(r["p_pred"]) > self.thresh
                and len(set(r["pfam"]) & bio_pfams) > self.n_biopfams)
                else 0
        }
        feature_cols = self.grouping + ["start", "end", target_col]
        feature_df = (pfam_df  .groupby(self.grouping, sort=False)
                            .apply(lambda x: x.assign(**cluster_condition))
                            .reset_index(drop=True)
                            .groupby(feature_cols, sort=False)
                            .first()
                            .reset_index()
                            .loc[:, feature_cols]
        )
        print(feature_df)
        clusters = [self._extract_clusters(df, target_col=target_col, prefix=s)
            for s, df in feature_df.groupby(self.seq_col)]

        print(clusters)

        return pfam_df

    def _extract_clusters(self, pfam_df, target_col="Y_pred", prefix="cluster"):
        cluster_num = 1
        cluster_state = False
        cluster_dict = {}
        for n in range(len(pfam_df)):
            row = pfam_df.iloc[n]
            if row[target_col] == 1:
                # non-cluster -> cluster
                if not cluster_state:
                    cluster_name = prefix + "_" + str(cluster_num)
                    row = (pd   .DataFrame(row[self.grouping + [n, "start", "end"]])
                                .transpose())
                    cluster_dict[cluster_name] = row
                    cluster_state = True
                # cluster -> cluster
                else:
                    cluster_dict[cluster_name] = cluster_dict[cluster_name].append(row)
            else:
                # cluster -> non-cluster
                if cluster_state:
                    cluster_num += 1
                    cluster_state = False
                # non-cluster -> non-cluster
                # do nothing
        return cluster_dict


    def _antismash_find_cf_clusters(self, pfam_df):
        """Find clusters based on ClusterFinder probabilities of CDS_motif features"""
        cf_clusters = []
        state = "seed"
        clusterpositions = [0,0]
        pfam_ids = []
        loop_index = 1
        probabilities = [0]
        pfam_df_arr = pfam_df.values
        features_df = pfam_df.groupby(self.grouping, sort=False).first().reset_index()
        for n, row in pfam_df.iterrows():
            featurepositions = int(row["start"]), int(row["end"])
            try:
                cf_probability = float(row["p_pred"])
            except:
                loop_index += 1
                continue
            if cf_probability >= 0.3:
                if state == "seed":
                    state = "extend"
                    probabilities = [cf_probability]
                    clusterpositions = [min(featurepositions), max(featurepositions)]
                    pfam_ids = []
                    pfam_ids.append(row["pfam"])
                    # pfam_ids.append([xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref][0].partition("PFAM: ")[2])
                else:
                    probabilities.append(cf_probability)
                    if max(featurepositions) > clusterpositions[1]:
                        clusterpositions[1] = max(featurepositions)
                    pfam_ids.append(row["pfam"])
                    # pfam_ids.append([xref for xref in feature.qualifiers['db_xref'] if "PFAM: " in xref][0].partition("PFAM: ")[2])
            else:
                if state == "extend":
                    state = "seed"
                    clusterpositions, cdsnr = self._antismash_find_nr_cds(row, features)
                    if self._antismash_is_good_cluster_hit(
                            cdsnr, probabilities, pfam_ids):

                        # logging.debug('Adding cluster %s at position %s to %s',
                        #     len(cf_clusters) + 1, clusterpositions[0],
                        #     clusterpositions[1])

                        cf_clusters.append((clusterpositions, mean(probabilities)))
                    clusterpositions = []
                    pfam_ids = []
            if loop_index == len(pfam_df):
                if len(clusterpositions) > 0:
                    clusterpositions, cdsnr = self._antismash_find_nr_cds(clusterpositions)
                    if self._antismash_is_good_cluster_hit(
                            cdsnr, probabilities, pfam_ids):

                        # logging.debug('Adding cluster %s at position %s to %s',
                        #     len(cf_clusters) + 1, clusterpositions[0],
                        #     clusterpositions[1])

                        cf_clusters.append((clusterpositions, mean(probabilities)))
                clusterpositions = []
                pfam_ids = []
            loop_index += 1
        return cf_clusters

    def _antismash_find_nr_cds(self, clusterpositions):
        """Find the number of CDSs in candidate cluster and adjust the cluster
        starts and ends to match the CDS starts and ends
        """
        cdsfeatures = utils.get_cds_features(seq_record)
        withinclustercdsfeatures = []
        for cds in cdsfeatures:
             if clusterpositions[0] <= int(cds.location.start) <= clusterpositions[1] or \
                clusterpositions[0] <= int(cds.location.end) <= clusterpositions[1] or \
                int(cds.location.start) <= clusterpositions[0] <= int(cds.location.end) or \
                int(cds.location.start) <= clusterpositions[1] <= int(cds.location.end):
                withinclustercdsfeatures.append(cds)
        if len(withinclustercdsfeatures) == 0:
            return clusterpositions, 0
        startlocations = [int(cds.location.start) for cds in withinclustercdsfeatures]
        endlocations = [int(cds.location.end) for cds in withinclustercdsfeatures]
        #If statement to avoid getting the complete genome as cluster if one CDS starts at end and finishes at start of genome
        if seq_record is not None and not (0 in startlocations and len(seq_record.seq) in endlocations):
            newclusterstart = min(startlocations)
            newclusterend = max(endlocations)
            newclusterpositions = [newclusterstart, newclusterend]
        else:
            newclusterpositions = clusterpositions
        return newclusterpositions, len(withinclustercdsfeatures)

    def _antismash_is_good_cluster_hit(self, cdsnr, probabilities, pfam_ids):
        """Check if the current cluster is a good hit"""

        biosynthetic_pfams = ["PF00109", "PF02801", "PF08659", "PF00378", "PF08541",
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
        "PF08484", "PF08421"]
        npfams = len(set(biosynthetic_pfams) & set(pfam_ids))
        rc = cdsnr >= options.cdsnr and \
             mean(probabilities) >= self.thresh and \
             npfams >= self.n_pfams