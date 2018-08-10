#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import time
import subprocess
import urllib.request
from setuptools import setup, find_packages

NAME = "ORION"
DESCRIPTION = "Predicting biosynthetic gene clusters using conditional random fields."
URL = "https://git.embl.de/fleck/ORION.git"
EMAIL = "jonas.simon.fleck@gmail.com"
AUTHOR = "Jonas Simon Fleck"
REQUIRES_PYTHON = ">=3.6.0"
VERSION = "0.1.0"
LICENSE = "GNU GPL 3.0"

# Python package requirements
REQUIRED = [
    "numpy >= 1.14.0",
    "scipy >= 1.1.0",
    "sklearn-crfsuite >= 0.3.6",
    "scikit-learn >= 0.19.1",
    "joblib >= 0.11",
    "pandas >= 0.22.0"
]

# Path to Pfam Database
PFAM_DB = "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz"

# Read README.md for long description
HERE = os.path.abspath(os.path.dirname(sys.argv[0]))

try:
    with open(os.path.join(HERE, "README.md"), "rt") as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# FUNC
def safe_extract(file):
    extract_cmd = "gzip -df " + file
    try:
        FNULL = open(os.devnull, 'w')
        process = subprocess.Popen(extract_cmd.split(),stderr=FNULL,stdout=FNULL)
        output, error = process.communicate()
    except:
        sys.stdout.write("Error: failed to extract files\n")
        sys.exit(1)
    if process.returncode:
        sys.stdout.write("Error: failed to extract files\n")
        sys.exit(1)
    else:
        sys.stdout.write("Done.\n")

def reporthook(count, block_size, total_size):
    """Reports download progress"""
    global start_time
    if count == 0:
        start_time = time.time()
        return
    duration = time.time() - start_time
    progress_size = int(count * block_size)
    speed = int(progress_size / (1024 * duration))
    percent = int(count * block_size * 100 / total_size)
    prog = int(50 * percent / 100)
    sys.stdout.write(
        "\r Pfam-A.hmm.gz: {0}{1} | {2}%, {3} MB, {4} KB/s, {5} seconds passed".format(
        "#" * prog, " " * (50 - prog),
        percent, int(progress_size / (1024 * 1024)), speed, int(duration)))
    sys.stdout.flush()


def query_yes_no(question, default="yes"):
    """Ask a yes/no question via raw_input() and return their answer.

    "question" is a string that is presented to the user.
    "default" is the presumed answer if the user just hits <Enter>.
        It must be "yes" (the default), "no" or None (meaning
        an answer is required of the user).

    The "answer" return value is True for "yes" or False for "no".
    """
    valid = {"yes": True, "y": True, "ye": True,
             "no": False, "n": False}
    if default is None:
        prompt = " (y/n) "
    elif default == "yes":
        prompt = " ([y]/n) "
    elif default == "no":
        prompt = " (y/[n]) "
    else:
        raise ValueError("invalid default answer: '{0}'".format(default))

    while True:
        sys.stdout.write(question + prompt)
        choice = input().lower()
        if default is not None and choice == '':
            return valid[default]
        elif choice in valid:
            return valid[choice]
        else:
            sys.stdout.write("Please respond with 'yes' or 'no' "
                             "(or 'y' or 'n').\n")

# MAIN
if __name__ == "__main__":

    sys.stdout.write(" ------------------------------------------------------\n")
    sys.stdout.write("|                     Greetings!                       |\n")
    sys.stdout.write("|            Welcome to the ORION installer            |\n")
    sys.stdout.write(" ------------------------------------------------------\n")

    default_dir = os.path.join(HERE, "data/pfam/")

    sys.stdout.write("Downloading the Pfam database (~1.2Gb).\n")
    default = query_yes_no("Do you want to download the database in the current program directory ({0})?".format(default_dir))

    if default:
        if not os.path.exists(default_dir):
            os.makedirs(default_dir)
        db_dir = default_dir
    else:
        db_dir = input("Please specify the directory where you wish to download the Pfam database:\n> ")
        if not os.path.exists(db_dir):
            os.makedirs(db_dir)

    db_path = os.path.join(os.path.realpath(db_dir), "Pfam-A.hmm.gz")

    force = False
    if os.path.exists(db_path[:-3]):
        force = query_yes_no("There is already a Pfam database in that directory. Do you want to overwrite it with the current download?", default="no")

    if force or not os.path.exists(db_path[:-3]):
        urllib.request.urlretrieve(PFAM_DB, db_path, reporthook)
        sys.stdout.write("\nDone.\n")
        sys.stdout.write("Decompressing the Pfam database.\n")
        safe_extract(db_path)

    config_path = os.path.join(HERE, "data/db_config.txt")
    with open(config_path, "wt") as f:
        f.write(db_path[:-3] + "\n")

    sys.stdout.write("Setting up Python dependencies\n\n")

    try:
        setup(
            name = NAME,
            version = VERSION,
            url = URL,
            author = AUTHOR,
            license = LICENSE,
            author_email = EMAIL,
            long_description = long_description,
            description = DESCRIPTION,
            python_requires = REQUIRES_PYTHON,
            packages = find_packages(),
            install_requires = REQUIRED,
        )
        sys.stdout.write("Done.\n")

    except:
        sys.stdout.write("It seems like the python dependencies could not be set up correctly. You might not have enough permissions to write to the PYTHONPATH.\nIf you are on a remote server, you can try re-running the installation from within a conda environment with pip installed.\nOtherwise, you might have to supply the following python dependencies manually:\n\n{0}\n\n".format(
        "\n".join(REQUIRED)))

    to_path = query_yes_no("Do you want to add orion.py to your PATH?")
    if to_path:
        os.system("chmod +x " + os.path.join(HERE, "orion.py"))
        os.system("echo '# ORION' >> ~/.bashrc")
        os.system("echo 'export PATH=\"{0}:$PATH\"' >> ~/.bashrc".format(HERE))
        os.system("source ~/.bashrc")

    sys.stdout.write("\nCongratulations! ORION has been successfully installed.\n\n")
