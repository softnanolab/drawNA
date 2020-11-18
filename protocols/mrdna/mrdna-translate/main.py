#!/usr/bin/env python
"""Example protocol where mrDNA is used to simulate replicas"""
from imports import *

import os

from softnanotools.runner import Runner
from softnanotools.logger import Logger
logger = Logger(__name__)

from backend import Manager

class Job(Runner):
    def __init__(self, manager: Manager):
        logger.info('Creating Runner')
        self.manager = manager
        return

    @Runner.task(0)
    def setup(self):
        logger.info('Running setup...')
        return

    @Runner.task(1)
    def loop(self):
        logger.info('Running loop...')
        return

    @Runner.task(2)
    def finalise(self):
        logger.info('Running finalise...')
        return

def main(**kwargs):
    
    if kwargs['tacoxDNA'] == None:
        logger.error(
            'tacoxDNA has not been set! It can be set by doing either:\n'
            '\t1. export TACOXDNA=/path/to/tacoxDNA\n'
            '\t2. python main.py ... --tacoxDNA /path/to/tacoxDNA'
        )
    manager = Manager(**kwargs)

    
    job = Job(manager)

    # Runs all methods decorated with @Runner.task() in order
    job.execute()
    return

if __name__ == "__main__":
    from argparse import ArgumentParser
    parser = ArgumentParser()
    parser.add_argument("-n", "--number", type=int, required=False, default=5, help='Number of replicas')
    parser.add_argument("-s", "--spacing", type=int, required=False, default=5, help='Spacing between replicas')
    parser.add_argument("--tacoxDNA", default=os.environ.get('TACOXDNA', None), help='Path to tacoxDNA root directory')
    main(**vars(parser.parse_args()))