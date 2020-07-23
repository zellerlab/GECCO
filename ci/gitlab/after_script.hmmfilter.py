import gzip
import re
import os
import io

import tqdm
from gecco.crf import ClusterCRF
from gecco.hmmer import embedded_hmms

# Create the artifact folder
os.makedirs(os.path.join("ci", "artifacts"), exist_ok=True)

# Load the internal CRF model and compile a regex that matches the domains
# known to the CRF (i.e. useful domains for us to annotate with)
crf = ClusterCRF.trained()
rx = re.compile("|".join(crf.model.attributes_).encode("utf-8"))

# Filter the hmms
for hmm in embedded_hmms():
    input_ = os.path.join("ci", "artifacts", "{}.hmm.gz".format(hmm.id))
    with gzip.open(hmm.path, "rb") as src, gzip.open(input_, "wb") as dst:
        blocklines = []
        for line in src:
            blocklines.append(line)
            if line == b"//\n":
                if any(rx.search(line) is not None for line in blocklines):
                    dst.writelines(blocklines)
                blocklines.clear()
