#!/usr/bin/env python
"""In this protocol, we generate single and double strands and show that the parameters for their generation are optimised."""

import subprocess
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from softnanotools.runner import Runner
from softnanotools.logger import Logger
logger = Logger(__name__)

from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.oxdna.strand import generate_helix

from quick_plot import main as _quick_plot # type: ignore

def get_params(ID: str) -> dict:
    params = {
        'conf_file': f'oxdna.{ID}.conf',
        'topology': f'oxdna.{ID}.top',
        'energy_file': f'oxdna.{ID}.energy',
        'trajectory_file': f'oxdna.{ID}.traj',
    }
    return params


class Benchmarks(Runner):
    def __init__(
        self, 
        oxDNA: str = 'oxDNA',
        input_file: str = 'oxdna.main.in',
        quick_plot: bool = False
    ):
        super().__init__()
        self.N = 21
        self.box = [30, 30, 30]
        self.oxDNA = oxDNA
        self.input_file = input_file
        self.quick_plot = quick_plot

    def run_simulation(
        self, 
        system: System, 
        params: dict,
    ) -> float:
        """Runs a generic oxDNA system and returns the difference in energy
        from the first step of the simulation and the average of the last 20.
        """
        if Path(params['energy_file']).exists():
            Path(params['energy_file']).unlink()
        system_name = params['conf_file'].split('.')[1]
        system.write_oxDNA(system_name)
        cmd = [self.oxDNA, self.input_file] 
        cmd += [f'{key}={value}' for key, value in params.items()]
        try:
            subprocess.check_output(cmd)
        except KeyboardInterrupt:
            logger.debug('Simulation interrupted by user using CTRL+C')
        return #self.evaluate_energy(params['energy_file'])

    def evaluate_energy(self, fname: str) -> float:
        """Returns the difference between the initial potential energy and
        final potential energy
        """
        return


    @Runner.task(0)
    def single_strand(self):
        """Runs a set of systems with 1 strand generated 
        using the drawNA package generate_helix function
        """
        system = System(self.box)
        strands = generate_helix(
            self.N, 
            sequence='A'*self.N
        )
        system.add_strands(strands)
        params = get_params('ss')
        self.run_simulation(system, params)
        if self.quick_plot:
            _quick_plot(params['energy_file'])
        return

    @Runner.task(1)
    def double_strand(self):
        """Runs a set of systems with 1 double strand generated 
        using the drawNA package generate_helix function
        """
        system = System(self.box)
        system = System(self.box)
        strands = generate_helix(
            self.N, 
            sequence='A'*self.N,
            double=True
        )
        system.add_strands(strands)
        params = get_params('ds')
        self.run_simulation(system, params)
        if self.quick_plot:
            _quick_plot(params['energy_file'])
        return

def main(**kwargs):
    Benchmarks(**kwargs).execute()
    if kwargs['quick_plot']:
        plt.show()
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()

    # oxDNA argument
    parser.add_argument(
        '-o', 
        '--oxDNA', 
        default='oxDNA', 
        help='Path to oxDNA executable'
    )

    # input_file argument
    parser.add_argument(
        '-i', 
        '--input-file', 
        default='oxdna.main.in', 
        help='Path to oxDNA input file template'
    )

    # input_file argument
    parser.add_argument(
        '-q', 
        '--quick-plot', 
        action='store_true', 
        help='Show energy plots'
    )

    # run main
    main(**vars(parser.parse_args()))