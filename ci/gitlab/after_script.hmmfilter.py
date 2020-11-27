import contextlib
import gzip
import re
import os
import io
import sys

import tqdm
import pkg_resources
import pyhmmer
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
    domains = [ domain.strip().decode() for domain in f ]

# Filter the hmms
for hmm_db in embedded_hmms():

    in_ = os.path.join("ci", "cache", "{}.{}.hmm".format(hmm_db.id, hmm_db.version))
    out = os.path.join("ci", "artifacts", "{}.hmm".format(hmm_db.id))
    size = sum(1 for e in interpro.entries if e.source_database.upper().startswith(hmm_db.id.upper()))

    with open(out, "wb") as dst:
        for hmm in tqdm.tqdm(pyhmmer.plan7.HMMFile(in_), desc=hmm_db.id, total=size):
            if hmm_db.relabel(hmm.accession.decode()) in domains:
                hmm.write(dst)
                nwritten += 1
