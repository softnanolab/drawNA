#!/usr/bin/env python
"""Generate parallel systems at different points in time"""
import subprocess
import shutil

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from drawNA.readers import OXDNAReader
from drawNA.oxdna.strand import generate_helix
from drawNA.oxdna import Nucleotide, Strand, System

def calculate_box_size(N: int, spacing: float) -> np.ndarray:
    """Determine the box size required"""
    box = np.array([5., 5., 5.]) * spacing * N
    return box

def generate_system(box: np.ndarray) -> System:
    """Generate a simple oxDNA system"""
    system = System(box)
    strand = generate_helix(40, double=False)[0]
    new = []
    for nt in strand.nucleotides[10:30]:
        new.append(nt.make_across())
    new = Strand(nucleotides=new[::-1])
    system.add_strands([strand, new])    
    return system

def run_simulation(oxDNA: str, input_file: str, **kwargs):
    """Run a simulation using oxDNA"""
    command = [oxDNA, input_file]
    for key, value in kwargs.items():
        command.append(f'{key}={value}')
    subprocess.check_output(command)
    return

def main(number: int, oxDNA: str, input_file: str):
    DELTA = 5.
    MAIN = 'main'
    TEMPLATE = MAIN + '.{}'
    N = number
    SIMULATIONS = range(1, N)

    main_system = generate_system(
        calculate_box_size(N, DELTA)
    )

    clone_system = main_system.copy()
    clone_name = TEMPLATE.format(0)
    clone_system.write_oxDNA(clone_name)

    for i in SIMULATIONS:
        new_name = TEMPLATE.format(i)

        run_simulation(
            oxDNA, 
            input_file,
            conf_file=f'oxdna.{clone_name}.conf',
            topology=f'oxdna.{clone_name}.top',
            lastconf_file=f'oxdna.{new_name}.conf',
        )

        # copy old topology to new one
        shutil.copy(f'oxdna.{clone_name}.top', f'oxdna.{new_name}.top')

        # use predefined new_* to define filenames and then import
        # this overwrites the previously defined clone_* variables
        clone_name = new_name
        clone_system = OXDNAReader([f'oxdna.{new_name}.conf', f'oxdna.{new_name}.top']).system

        # do things with the clone system
        clone_system.translate(np.array([0., 0., DELTA]))

        # clone the strands and add them to the main system
        clone_strands = clone_system.strands
        main_system.add_strands(clone_strands)
        
        # standard operations  to export everything correctly
        # and run simulation
        clone_system.write_oxDNA(clone_name)

    main_system.write_oxDNA(MAIN)
    run_simulation(
        oxDNA, 
        input_file, 
        conf_file=f'oxdna.{MAIN}.conf',
        topology=f'oxdna.{MAIN}.top',
        steps=1000000,
    )
    shutil.copy(
        f'oxdna.{MAIN}.top', 
        f'oxdna.traj.top'
    )
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-n", "--number", type=int, required=False, default=5)
    parser.add_argument("-o", "--oxDNA", type=str, required=False, default='oxDNA')
    parser.add_argument("-i", "--input-file", type=str, required=False, default='input')
    main(**vars(parser.parse_args()))