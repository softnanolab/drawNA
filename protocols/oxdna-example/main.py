#!/usr/bin/env python
"""Update description for oxdna-example"""

def main(dummy : None = None):
    print("Automatically Generated Program for oxdna-example")
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("--dummy", type=None, required=False, default=None)
    main(**vars(parser.parse_args()))