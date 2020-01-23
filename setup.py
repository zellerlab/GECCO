#!/usr/bin/env python
# -*- coding: utf-8 -*-

import configparser
import hashlib
import io
import os
import sys
import urllib.request

import setuptools
from setuptools.command.build_ext import build_ext as _build_ext
from setuptools.command.sdist import sdist as _sdist


class Resource(setuptools.Extension):
    """A phony `Extension` that will download resource when "compiled".
    """

    def __init__(self, url, name):
        self.url = url
        self.target = os.path.join(*name.split('.'))
        setuptools.Extension.__init__(self, name, [f"{self.target}.hmm.gz.md5"])


class build_ext(_build_ext):
    """An extension to the `build_ext` command that can download resources.
    """

    def _download(self, src_url, dst_path):
        import tqdm

        hasher = hashlib.md5()
        with urllib.request.urlopen(src_url) as src:
            with open(dst_path, "wb") as dst:

                content_length = int(src.headers['Content-Length'])
                read = lambda: src.read(io.DEFAULT_BUFFER_SIZE)

                pbar = tqdm.tqdm(
                    total=content_length,
                    unit='B',
                    unit_scale=True,
                    unit_divisor=1024,
                    leave=False,
                    file=sys.stdout,
                    miniters=-1,
                )

                for chunk in iter(read, b''):
                    hasher.update(chunk)
                    pbar.update(dst.write(chunk))

                    if not sys.stdout.isatty():
                        pbar.refresh()
                        sys.stdout.write('\n')

        pbar.close()
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
        actual_md5 = self._download(ext.url, dst_file)

        #
        self.announce(f"checking MD5 checksum of {dst_file}")
        with open(ext.sources[0]) as f:
            expected_md5 = f.readline().strip()

        #
        assert actual_md5 == expected_md5


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


if __name__ == "__main__":
    setuptools.setup(
        cmdclass={
            "build_ext": build_ext,
            "sdist": sdist,
        },
        ext_modules=[
            Resource(
                "ftp://ftp.ebi.ac.uk/pub/databases/Pfam/releases/Pfam31.0/Pfam-A.hmm.gz",
                "gecco.data.hmms.Pfam-A"
            )],
    )
