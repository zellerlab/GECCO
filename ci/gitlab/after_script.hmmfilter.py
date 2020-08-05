import contextlib
import gzip
import re
import os
import io
import sys

import tqdm
import pkg_resources
sys.path.insert(0, os.path.realpath(os.path.join(__file__, "..", "..", "..")))

from gecco.hmmer import embedded_hmms
from gecco.interpro import InterPro

# Create the artifact folder
os.makedirs(os.path.join("ci", "artifacts"), exist_ok=True)

# Load InterPro to know how many entries we have to process
interpro = InterPro.load()

# Load the domains used by the CRF and compile a regex that matches the domains
# known to the CRF (i.e. useful domains for us to annotate with)
with pkg_resources.resource_stream("gecco.types", "domains.tsv") as f:
    domains = [ domain.strip() for domain in f ]
    rx = re.compile(b"|".join(domains))

# Filter the hmms
for hmm in embedded_hmms():
    out = os.path.join("ci", "artifacts", "{}.hmm.gz".format(hmm.id))
    in_ = os.path.join("ci", "cache", "{}.{}.hmm.gz".format(hmm.id, hmm.version))
    size = sum(1 for e in interpro.entries if e.source_database.upper().startswith(hmm.id.upper()))
    pbar = tqdm.tqdm(desc=hmm.id, total=size)

    with contextlib.ExitStack() as ctx:
        pbar = ctx.enter_context(pbar)
        src = ctx.enter_context(gzip.open(in_, "rb"))
        dst = ctx.enter_context(gzip.open(out, "wb"))

        blocklines = []
        for line in src:
            blocklines.append(line)
            if line == b"//\n":
                if any(rx.search(line) is not None for line in blocklines):
                    dst.writelines(blocklines)
                blocklines.clear()
                pbar.update(1)
