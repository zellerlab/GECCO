import argparse

# FUNC
def main_interface():
    parser = argparse.ArgumentParser(description="Predicts biosynthetic gene clusters from a genome FASTA file.")

    parser.add_argument("FASTA",
                        type=str,
                        metavar="<FASTA>",
                        help="A genome FASTA file as input.")

    parser.add_argument("-o", "--output-dir",
                        dest="out",
                        type=str,
                        default="./",
                        metavar="<out_dir>",
                        help="Output directory.")

    parser.add_argument("--e-filter",
                        dest="e_filter",
                        type=float,
                        default="1e-5",
                        metavar="<e_filter>",
                        help="E-value cutoff for pfam domains to be included.")

    parser.add_argument("--feature-type",
                        dest="feature_type",
                        default="protein",
                        type=str,
                        metavar="<feature_type>",
                        help="How features should be extracted. 'Single', 'overlap' or on some grouping level ('group').")

    parser.add_argument("--weight-type",
                        dest="weight_type",
                        type=str,
                        metavar="<weight_type>",
                        help="Type of local weights for features.")

    args = parser.parse_args()
    return args

def scripts_interface():
    parser = argparse.ArgumentParser(description="Generic interface for all scripts which use the ClusterCRF (orion_[cv/loto/train/predict].py)")

    parser.add_argument("DATA",
                        type=str,
                        metavar="<DATA>",
                        help="Pfam table.")

    parser.add_argument("-o", "--output-basename",
                        dest="out",
                        type=str,
                        default="CRF",
                        metavar="<basename>",
                        help="Basename for the output.")

    parser.add_argument("-m", "--model",
                        dest="model",
                        type=str,
                        metavar="<model>",
                        help="Model for predictions.")

    parser.add_argument("-t", "--threads",
                        dest="threads",
                        type=int,
                        metavar="<threads>",
                        help="Number of CPUs.")

    parser.add_argument("-e", "--evalue",
                        dest="e_filter",
                        type=float,
                        default="1e-5",
                        help="E-value threshold for the test set.")

    parser.add_argument("-y", "--y-col",
                        dest="y",
                        type=str,
                        default="BGC",
                        help="Column with class labels.")

    parser.add_argument("-w", "--weight-col",
                        dest="w",
                        type=str,
                        nargs="+",
                        default=["1"],
                        help="Column to be used as local weights on features.")

    parser.add_argument("-f", "--feature-col",
                        dest="feat",
                        type=str,
                        nargs="+",
                        default=["pfam"],
                        help="Column to be used as features.")

    parser.add_argument("-s", "--split-col",
                        dest="split_col",
                        default="genome_id",
                        type=str,
                        help="Column to be used for splitting in to samples, i.e. different sequences.")

    parser.add_argument("-g", "--group-col",
                        dest="group_col",
                        default="protein_id",
                        type=str,
                        help="Column to be used for grouping features if feature_type is 'group'.")

    parser.add_argument("-p", "--threshold",
                        dest="thresh",
                        default="0.6",
                        type=float,
                        help="Probability threshold for clusters prediction.")

    parser.add_argument("--sort-cols",
                        dest="sort_cols",
                        default=["genome_id", "start", "domain_start"],
                        nargs="+",
                        type=str,
                        help="Columns to be used for sorting the data.")

    parser.add_argument("--strat-col",
                        dest="strat_col",
                        default="BGC_type",
                        type=str,
                        help="Columns to be used for stratifying the samples (BGC types).")

    parser.add_argument("--feature-type",
                        dest="feature_type",
                        type=str,
                        default="single",
                        help="How features should be extracted. 'Single', 'overlap' or on some grouping level ('group').")

    parser.add_argument("--truncate",
                        dest="truncate",
                        type=int,
                        help="Training set will be truncated to this length.")

    parser.add_argument("--overlap",
                        dest="overlap",
                        type=int,
                        default="2",
                        help="If overlapping features: How much overlap.")

    parser.add_argument("--no-shuffle",
                        dest="shuffle",
                        action="store_false",
                        help="Switch to turn of shuffling of the data before doing CV.")

    parser.add_argument("--folds",
                        dest="splits",
                        default="10",
                        type=int,
                        help="Number of folds for CV.")

    parser.add_argument("--C1",
                        dest="C1",
                        default="0.15",
                        type=float,
                        help="Parameter for L1 regularization.")

    parser.add_argument("--C2",
                        dest="C2",
                        default="1.75",
                        type=float,
                        help="Parameter for L2 regularization.")

    args = parser.parse_args()
    return args
