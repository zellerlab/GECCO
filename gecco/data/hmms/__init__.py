import configparser
import glob
import re
import typing
import pkg_resources

class Hmm(typing.NamedTuple):

    id: str
    version: str
    url: str
    action: str
    md5: typing.Optional[str] = None
    matching: typing.Optional[str] = None
    relabel_with: typing.Optional[str] = None

    @property
    def path(self):
        basename = f"{self.id}.hmm.gz"
        return pkg_resources.resource_filename(__name__, basename)

    def relabel(self, domains):
        if self.relabel_with is None:
            return domains
        before, after = re.match("^s/(.*)/(.*)/$", self.relabel_with).groups()
        regex = re.compile(before)
        return [regex.sub(after, domain) for domain in domains]


def iter():
    for ini in glob.glob(pkg_resources.resource_filename(__name__, "*.ini")):
        cfg = configparser.ConfigParser()
        cfg.read(ini)
        yield Hmm(**dict(cfg.items("hmm")))
