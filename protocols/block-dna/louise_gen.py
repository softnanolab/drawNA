#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 17:43:05 2020

@author: louiserosset
"""
from argparse import ArgumentParser
import numpy as np

from fiona_gen import generate_system
from system import SystemWithFolder

if __name__ == '__main__':
    parser = ArgumentParser()
    parser.add_argument('-n', '--number', type=int, nargs=3, required=True)
    parser.add_argument('-ds', '--double-stranded', type=int, nargs=3, required=True)
    parser.add_argument('-ss', '--single-stranded', type=int, nargs=3,required=True)
    parser.add_argument('-f', '--output-prefix')
    args = parser.parse_args()
    
    tot_length = np.arange(args.number[0], args.number[2]+1, args.number[1])
    ds_length = np.arange(args.double_stranded[0], args.double_stranded[2]+1, args.double_stranded[1])
    ss_length = np.arange(args.single_stranded[0], args.single_stranded[2]+1, args.single_stranded[1])

    sim_nb = 0  

    for i in range (len(tot_length)):
        for j in range (len(ds_length)):
            for k in range (len(ss_length)):
                system = generate_system(tot_length[i], ds_length[j], ss_length[k])
                sim_nb +=1

                if args.output_prefix:
                    system.write_oxDNA_folder(sim_nb, tot_length[i], ds_length[j], ss_length [k], args.output_prefix)