#!/usr/bin/env python
"""Example protocol where mrDNA is used to simulate replicas"""

from imports import *

class Job(softnanotools.runner.Runner):
    def __init__(self):
        return

def main(dummy : None = None):
    print("Automatically Generated Program for mrdna-translate")
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--dummy", type=None, required=False, default=None)
    main(**vars(parser.parse_args()))