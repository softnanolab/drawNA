#!/usr/bin/env python

import pandas as pd
import matplotlib.pyplot as plt

def main(energy: str, out: str = '', show: bool = False):
    fig, ax = plt.subplots()
    ax.set_xlabel('time')
    ax.set_ylabel('energy')
    data = pd.read_csv(energy, delim_whitespace=True, header=None)
    moving = data.rolling(15, center=True).mean()
    ax.plot(data[0], data[1])
    ax.plot(moving[0], moving[1])
    ax.set_xscale('log')
    ax.set_title(energy)
    if out:
        fig.savefig(out)
    if show:
        plt.show()
    return

if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('energy')
    parser.add_argument('-o', '--out', default='')
    parser.add_argument('--show', action='store_true')
    main(**vars(parser.parse_args()))

