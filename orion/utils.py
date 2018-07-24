

def convert_hmmer(dom_file, out_file):
    """Converts HMMER --domtblout output to regular TSV"""

    header = ["protein_id", "pfam", "i_Evalue", "domain_start", "domain_end"]

    with open(dom_file, "rt") as f:
        with open(out_file, "wt") as fout:
            fout.write("\t".join(header) + "\n")
            for line in f:
                if not line.startswith("#"):
                    line = line.split()
                    row = [line[0]] + [line[4]] + [line[12]] + line[17:19]
                    fout.write("\t".join(row) + "\n")
