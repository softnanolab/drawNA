#!/usr/bin/env python
"""Example protocol where mrDNA is used to simulate replicas"""
from imports import *

import os
import shutil
import numpy as np

from pathlib import Path

from softnanotools.runner import Runner
from softnanotools.logger import Logger
logger = Logger(__name__)

from backend import Manager

class Job(Manager, Runner):
    def __init__(self, **kwargs):
        logger.info('Creating Runner...')
        Manager.__init__(self, **kwargs)
        logger.info('Success!')
    
    @Runner.task(0)
    def initialise(self):
        logger.info('Running Initialisation...')
        self.define_multiple_systems()
        self.generate_system()
        self.generate_forces()
        self.generate_input_file('equilibration')
        logger.info('Success!')
        return

    @Runner.task(1)
    def equilibrate(self):
        logger.info('Running Equilibration...')
        self.run_equilibration()
        logger.info('Success!')
        return

    @Runner.task(2)
    def replicate(self):
        logger.info('Running Replication...')
        # run simulation ensuring that
        # that the input file generated has
        # the correct settings to allow
        # for the correct number of replicas
        self.generate_input_file('replication')
        self.run_replication()
        self.split_trajectory()
        logger.info('Success!')
        return

    @Runner.task(3)
    def compile(self):
        logger.info('Running Compilation...')
        # read trajectory and separate into
        # different files for each replica

        # create a new empty system
        # each replica in a new lattice point
        # in the system

        self.create_mother_system()

        # write the system to oxDNA and convert to PDB
        logger.info('Success!')
        return

    @Runner.task(4)
    def run(self):
        logger.info('Running Master Simulation...')
        # run the mrdna simulation using the previously
        # created PDB
        self.mrdna_simulate()
        # tidy up any files
        #self.cleanup_files()
        logger.info('Success!')
        return

    @Runner.task(5)
    def move(self):
        logger.info('Moving everything...')
        self.move_sim_files()
        logger.info('Success!')
        return






def main(args):
    total_length = np.arange(args['overall_length'][0], args['overall_length'][2]+1, args['overall_length'][1])
    percentage = np.arange(args['percent_stapled'][0], args['percent_stapled'][2]+1, args['percent_stapled'][1])
    stapled = np.arange(args['length_stapled'][0], args['length_stapled'][2]+1, args['length_stapled'][1])

    sim_index = -1

    for this_length in total_length:
        for this_percentage in percentage:
            for this_stapled in stapled:
                kwargs = args.copy()
                kwargs['overall_length'] = this_length
                kwargs['percent_stapled'] = this_percentage
                kwargs['length_stapled'] = this_stapled
                kwargs['box'] = args['box']

                sim_index += 1
                kwargs['sim_index'] = sim_index
          
                # setup Job
                job = Job(**kwargs)

                # execute Job
                job.execute()
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-b", "--box", type=float, default=40., help='Length of cubic box')
    parser.add_argument("-n", "--number", type=int, default=8, help='Number of replicas')
    #parser.add_argument("-ns", "--n_strands", type=int, default=8, help='Number of strands within one strand (defines number of turns)')
    parser.add_argument("-l", "--overall_length", type=int, nargs=3, default=8, help='Number of strands within one strand (defines number of turns)')
    #parser.add_argument("-l", "--length", type=int, default=16, help='Length of a single-stranded portion in between turns')
    parser.add_argument("-p", "--percent_stapled", type=int, nargs=3, default=16, help='Length of a single-stranded portion in between turns')
    #parser.add_argument("-s", "--stapled", type=int, default=5, help='Length of double-stranded portion')
    parser.add_argument("-ds", "--length_stapled", type=int, nargs=3, default=5, help='Length of double-stranded portion')
    parser.add_argument("--oxDNA", default='oxDNA', help='Path to oxDNA binary executable')
    parser.add_argument("--tacoxDNA", default=os.environ.get('TACOXDNA', None), help='Path to tacoxDNA root directory')
    parser.add_argument("--name", default='main', help='Naming-pattern for all of the resulting files')
    parser.add_argument("--root", default='data', help='Results Directory (will be created if it does not exist)')
    main(vars(parser.parse_args()))