import os
import sys
import warnings
warnings.filterwarnings("ignore", message="numpy.dtype size changed")
warnings.filterwarnings("ignore", message="numpy.ufunc size changed")
GECCO = os.path.abspath(os.path.dirname(os.path.abspath(sys.argv[0])) + "/..")
sys.path.append(GECCO)
import numpy as np
import pandas as pd
from gecco.hmmer import HMMER
from gecco.orf import ORFFinder
from gecco.interface import annot_interface

### TEST ###
# python gecco_annotate.py -p ../test/test.faa --db /g/scb2/zeller/fleck/DB/Pfam-A.hmm  -o /g/scb2/zeller/fleck/test/


# MAIN
if __name__ == "__main__":
    # PARAMS
    args = annot_interface()

    # Make out directory
    out_dir = args.out
    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    e_filter = min(1, args.e_filter)

    if args.GENOME:
        genome = args.GENOME
        base = ".".join(os.path.basename(genome).split(".")[:-1])

        # PRODIGAL
        prodigal_out = os.path.join(out_dir, "prodigal/")
        if not os.path.exists(prodigal_out):
            os.makedirs(prodigal_out)

        # Extract ORFs from genome
        prodigal = ORFFinder(genome, prodigal_out, method="prodigal")
        orf_file = prodigal.run()
        prodigal = True

    else:
        orf_file = args.PROTEINS
        base = ".".join(os.path.basename(orf_file).split(".")[:-1])
        prodigal = False

    # HMMER
    hmmer_out = os.path.join(out_dir, "hmmer/")
    if not os.path.exists(hmmer_out):
        os.makedirs(hmmer_out)

    # Run PFAM HMM DB over ORFs to annotate with HMM DB
    hmmer = HMMER(orf_file, hmmer_out, hmms=args.DB, prodigal=prodigal)
    pfam_df = hmmer.run()

    # Filter i-Evalue
    pfam_df = pfam_df[pfam_df["i_Evalue"] < e_filter]
    # Reformat pfam IDs
    pfam_df = pfam_df.assign(
        pfam = pfam_df["pfam"].str.replace(r"(PF\d+)\.\d+", lambda m: m.group(1)),
        sequence_id = base
    )

    # Write feature table to file
    feat_out = os.path.join(out_dir, base + ".features.tsv")
    pfam_df.to_csv(feat_out, sep="\t", index=False)
