import contextlib
import gzip
import re
import os
import io

import tqdm
from gecco.crf import ClusterCRF
from gecco.hmmer import embedded_hmms
from gecco.interpro import InterPro

# Create the artifact folder
os.makedirs(os.path.join("ci", "artifacts"), exist_ok=True)

# Load InterPro to know how many entries we have to process
interpro = InterPro.load()

# Load the internal CRF model and compile a regex that matches the domains
# known to the CRF (i.e. useful domains for us to annotate with)
crf = ClusterCRF.trained()
rx = re.compile("|".join(crf.model.attributes_).encode("utf-8"))


# Filter the hmms
for hmm in embedded_hmms():
    out = os.path.join("ci", "artifacts", "{}.hmm.gz".format(hmm.id))
    size = sum(1 for e in interpro.entries if e.source_database.upper().startswith(hmm.id.upper()))
    pbar = tqdm.tqdm(desc=hmm.id, total=size, unit="B", unit_scale=True, unit_divisor=1024)

    with contextlib.ExitStack() as ctx:
        pbar = ctx.enter_context(pbar)
        src = ctx.enter_context(gzip.open(hmm.path, "rb"))
        dst = ctx.enter_context(gzip.open(out, "wb"))

        blocklines = []
        for line in src:
            blocklines.append(line)
            if line == b"//\n":
                if any(rx.search(line) is not None for line in blocklines):
                    dst.writelines(blocklines)
                blocklines.clear()
                pbar.update(1)
