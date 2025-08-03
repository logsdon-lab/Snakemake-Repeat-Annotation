import os
import sys
import yaml
import argparse


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--infiles", nargs="+")
    ap.add_argument("-c", "--config", type=argparse.FileType("rb"), required=True)
    args = ap.parse_args()

    cfg = yaml.safe_load(args.config)

    cfg["samples"] = {}

    for file in args.infiles:
        bname, _ = os.path.splitext(os.path.basename(file).replace(".gz", ""))
        cfg["samples"][bname] = {"fa": file}

    yaml.safe_dump(cfg, sys.stdout)


if __name__ == "__main__":
    raise SystemExit(main())
