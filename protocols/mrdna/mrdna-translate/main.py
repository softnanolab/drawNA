#!/usr/bin/env python
"""Example protocol where mrDNA is used to simulate replicas"""
from imports import *

import os
import shutil

from pathlib import Path

from softnanotools.runner import Runner
from softnanotools.logger import Logger
logger = Logger(__name__)

from backend import Manager

class Job(Runner, Manager):
    def __init__(self, **kwargs):
        logger.info('Creating Runner')
        Runner.__init__(self)
        Manager.__init__(self, **kwargs)

    @Runner.task(0)
    def setup(self):
        """Setup single instance of a system"""
        logger.info('Running setup...')

        # create system
        self.system = generate_system()

        # create spherical force
        self.initialise_forces()
        return

    @Runner.task(1)
    def equilibrate(self):
        """Equilibrate single oxDNA system"""
        logger.info('Running equilibration...')
        equilibrated = False

        # run one simulation for a large number of steps
        self.oxdna_simulation()
        while not equilibrated:

            # run simulation for a shorter time
            self.oxdna_simulate(self.input_file)

            # check energy
            if self.check_energy():
                equilibrated = True

        return

    @Runner.task(2)
    def replicate(self):
        """Run a simulation that outputs a different """
        logger.info('Running replication...')
        return

    @Runner.task(3)
    def assemble(self):
        """Iterate through all simulations and create a large, compiled one"""
        logger.info('Running assemble')
        return

    @Runner.task(4)
    def finalise(self):
        """Create mrdna simulation and run to completion"""
        return

def main(**kwargs):
    
    # check tacoxDNA is installed properly
    if kwargs['tacoxDNA'] == None:
        logger.error(
            'tacoxDNA has not been set! It can be set by doing either:\n'
            '\t1. export TACOXDNA=/path/to/tacoxDNA\n'
            '\t2. python main.py ... --tacoxDNA /path/to/tacoxDNA'
        )
    
    # setup manager and job
    job = Job(**kwargs)

    # Runs all methods decorated with @Runner.task() in order
    job.execute()
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-n", "--number", type=int, required=False, default=5, help='Number of replicas')
    parser.add_argument("-s", "--spacing", type=int, required=False, default=5, help='Spacing between replicas')
    parser.add_argument("--tacoxDNA", default=os.environ.get('TACOXDNA', None), help='Path to tacoxDNA root directory')
    parser.add_argument("--oxDNA", default='oxDNA', help='Path to tacoxDNA root directory')
    parser.add_argument("-i", "--input-file", default='input', help='Path to tacoxDNA root directory')
    main(**vars(parser.parse_args()))