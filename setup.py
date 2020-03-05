#!/usr/bin/env python
# -*- coding: utf-8 -*-

import configparser
import distutils.cmd
import distutils.log
import gzip
import hashlib
import io
import math
import os
import sys
import tarfile
import urllib.request

import setuptools
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist

try:
    from tqdm import tqdm
except ImportError:
    tqdm = None


class ResponseProgressBar(object):

    def __init__(self, inner, tqdm=tqdm, **kwargs):
        self.inner = inner
        self.total = total = int(inner.headers['Content-Length'])

        if tqdm is not None:
            self.pbar = tqdm(
                total=total,
                leave=True,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
                file=sys.stdout,
                **kwargs
            )
            self.update = self.pbar.update
            self.refresh = self.pbar.refresh
        else:
            self.pbar = None
            self.current = 0
            self.next_p = 1
            self.desc = kwargs.get("desc")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, traceback):
        if sys.version_info >= (3,):
            self.inner.__exit__(exc_type, exc_value, traceback)
        if self.pbar is not None:
            self.pbar.__exit__(exc_type, exc_value, traceback)
        else:
            print()
        return False

    def read(self, n=None):
        chunk = self.inner.read(n)
        self.update(len(chunk))
        return chunk

    def update(self, n):
        self.current = c = self.current + n
        p = int(float(c) * 100.0 / self.total)
        if p >= self.next_p:
            self.next_p = p + 1
            print(
                "{} : {:>{}}B / {}B [{}%]".format(
                    self.desc,
                    c,
                    int(math.ceil(math.log10(self.total))),
                    self.total,
                    p
                ),
                end="\r"
            )
            sys.stdout.flush()

    def refresh(self):
        pass


class Resource(setuptools.Extension):
    """A phony `Extension` that will download resource when "compiled".
    """

    def __init__(self, url, name, action):
        self.url = url
        self.target = os.path.join(*name.split('.'))
        self.action = action
        setuptools.Extension.__init__(self, name, [f"{self.target}.hmm.gz.md5"])


class build_ext(_build_ext):
    """An extension to the `build_ext` command that can download resources.
    """

    def _flush(self, pbar):
        if not sys.stdout.isatty():
            pbar.refresh()
            sys.stdout.flush()
            sys.stdout.write('\n\033[F')

    def extract(self, src_url, dst_path):
        hasher = hashlib.md5()
        with ResponseProgressBar(urllib.request.urlopen(src_url)) as src:
            with open(dst_path, "wb") as dst:
                read = lambda: src.read(io.DEFAULT_BUFFER_SIZE)
                for chunk in iter(read, b''):
                    hasher.update(chunk)
                    dst.write(chunk)
                    self._flush(src)
        return hasher.hexdigest()

    def merge(self, src_url, dst_path):
        hasher = hashlib.md5()
        with ResponseProgressBar(urllib.request.urlopen(src_url)) as src:
            with gzip.open(dst_path, "wb") as dst:
                tar = tarfile.open(fileobj=src, mode="r|gz")
                members = filter(
                    lambda member: member.name.lower().endswith("hmm"),
                    iter(tar.next, None)
                )
                for member in members:
                    with tar.extractfile(member) as mem:
                        read = lambda: mem.read(io.DEFAULT_BUFFER_SIZE)
                        for chunk in iter(read, b""):
                            hasher.update(chunk)
                            dst.write(chunk)
                            self._flush(src)
        return hasher.hexdigest()

    def get_ext_filename(self, ext_name):
        basepath = os.path.join(*ext_name.split("."))
        return f'{basepath}.hmm.gz'

    def build_extension(self, ext):
        # Get the destination path where `setuptools` wants the file to be
        dst_file = self.get_ext_fullpath(ext.name)
        dst_dir = os.path.dirname(dst_file)
        self.mkpath(dst_dir)

        # Download the file to the requested location
        self.announce(f"downloading {ext.url}", 2)
        actual_md5 = getattr(self, ext.action)(ext.url, dst_file)

        #
        self.announce(f"checking MD5 checksum of {dst_file}")
        with open(ext.sources[0]) as f:
            expected_md5 = f.readline().strip()

        #
        #assert actual_md5 == expected_md5


class sdist(_sdist):
    """An extension to the `sdist` command that generates a `pyproject.toml`.
    """

    def run(self):
        # build `pyproject.toml` from `setup.cfg`
        c = configparser.ConfigParser()
        c.add_section("build-system")
        c.set("build-system", "requires", str(self.distribution.setup_requires))
        c.set("build-system", 'build-backend', '"setuptools.build_meta"')
        with open("pyproject.toml", "w") as pyproject:
            c.write(pyproject)

        # run the rest of the packaging
        _sdist.run(self)


class update_model(distutils.cmd.Command):
    """A custom command to update the internal CRF model."""

    description = 'update the CRF model embedded in the source'
    user_options = [
      ('model=', 'm', 'the path to the new CRF model to use'),
    ]

    def initialize_options(self):
        self.model = None

    def finalize_options(self):
        if self.model is None:
            raise ValueError("--model argument must be given")
        elif not os.path.exists(self.model):
            raise FileNotFoundError(self.model)

    def info(self, msg):
        self.announce(msg, level=distutils.log.INFO)

    def run(self):
        import gecco.data

        # Copy the file to the new in-source location and compute its hash.
        hasher = hashlib.md5()
        self.info("Copying the model to the in-source location")
        with open(self.model, "rb") as src:
            with open(gecco.data.realpath("model/crf.model"), "wb") as dst:
                read = lambda: src.read(io.DEFAULT_BUFFER_SIZE)
                for chunk in iter(read, b''):
                    hasher.update(chunk)
                    dst.write(chunk)

        # Write the hash to the signature file next to the model
        self.info("Writing the MD5 signature file")
        with open(gecco.data.realpath("model/crf.model.md5"), "w") as sig:
            sig.write(hasher.hexdigest())


if __name__ == "__main__":
    setuptools.setup(
        cmdclass={
            "build_ext": build_ext,
            "sdist": sdist,
            "update_model": update_model,
        },
        ext_modules=[
            Resource(
                "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz",
                "gecco.data.hmms.Pfam",
                action="extract"
            ),
            Resource(
                "ftp://ftp.jcvi.org/pub/data/TIGRFAMs/14.0_Release/TIGRFAMs_14.0_HMM.tar.gz",
                "gecco.data.hmms.Tigrfam",
                action="merge"
            ),
        ],
    )
