import argparse
import configparser


parser = argparse.ArgumentParser()
parser.add_argument("-i", "--input", required=True)
parser.add_argument("-o", "--output", required=True)
args = parser.parse_args()

cfg = configparser.ConfigParser()
cfg.read(args.input)

with open(args.output, "w") as f:
    f.write(cfg.get("options", "setup_requires"))
    f.write(cfg.get("options", "install_requires"))
    for _, v in cfg.items("options.extras_require"):
        f.write(v)