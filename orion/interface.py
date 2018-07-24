import argparse

# FUNC
def refine_interface():
    parser = argparse.ArgumentParser(description="Interface for the refinement script.")

    parser.add_argument("DATA",
                        type=str,
                        metavar="<DATA>",
                        help="Pfam table.")

    parser.add_argument("-o", "--output-basename",
                        dest="out",
                        type=str,
                        default="CRF",
                        metavar="<basename>",
                        help="Basename for predictions.")

    parser.add_argument("-t", "--threshold",
                        dest="thresh",
                        default="0.5",
                        type=float,
                        help="Probability threshold for clusters prediction.")

    args = parser.parse_args()
    return args

def crf_interface():
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
                        help="Basename for the output table.")

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
                        default="1",
                        help="E-value threshold for the test set.")

    parser.add_argument("-w", "--weight-col",
                        dest="w",
                        default="rev_i_Evalue",
                        type=str,
                        help="Column to be used as local weights on pfam domains.")

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

    args = parser.parse_args()
    return args
