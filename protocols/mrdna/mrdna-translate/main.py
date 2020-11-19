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

class Job(Runner):
    def __init__(self, manager: Manager):
        logger.info('Creating Runner')
        self.manager = manager
        return

    @Runner.task(0)
    def setup(self):
        logger.info('Running setup...')
        system = self.manager.generate_system()
        self.manager.cache['pdb'] = self.manager.export_PDB(
            system, 
            f'{self.manager.name}.0'
        )
        self.manager.child_systems[0] = system
        # run mrDNA simulation [0]
        self.manager.simulate(
            self.manager.cache['pdb']
        )
        return

    @Runner.task(1)
    def loop(self):
        logger.info('Running loop...')
        try:
            for i, system in enumerate(self.manager.child_systems[1:]):
                # import resulting PDB file as 
                # oxDNA system [1]
                system = self.manager.import_PDB(self.manager.cache['pdb'])

                # recentre [1]
                system = self.manager.recenter_system(system)

                # export PDB file [1]
                self.manager.cache['pdb'] = self.manager.export_PDB(
                    system,
                    f'{self.manager.name}.{i}'
                )
                logger.debug(f"cache['pdb']: {self.manager.cache['pdb']}")
                assert Path(self.manager.cache['pdb']).exists()
                # run simulation [1]
                self.manager.simulate(
                    str(Path(self.manager.cache['pdb']).resolve())
                )
        except:
            logger.warning('Loop Failed')
        return

    @Runner.task(2)
    def finalise(self):
        logger.info('Running finalise...')
        # for each child system
        for i, system in enumerate(self.manager.child_systems):
            # translate a copy of the system
            temp_system = self.manager.recenter_system(system.copy())
            # export child system
            temp_system.write_oxDNA(f'{self.manager.name}.{i}', root=self.manager.root)
            temp_system.translate(self.manager.spacing)

            # take the strands and
            # add them to the mother system
            for strand in temp_system.strands:
                self.manager.mother_system.add_strand(
                    strand.copy()
                )

        # export mother system
        self.manager.mother_system.write_oxDNA(self.manager.name, root=self.manager.root)

        # delete other files
        for ext in [
            "*.txt*", 
            "*.bd*", 
            "*.pdb*", 
            "*.psf*", 
            "*.namd*",
            "oxdna*",
            "*.exb",
        ]:
            logger.info(f'Deleting {ext}')
            for p in Path(".").glob(ext):
                p.unlink()

        for folder in ['charmm36.nbfix', 'output', 'potentials']:
            logger.info(f'Deleting {folder}')
            subprocess.check_output(['rm', '-rf', folder])

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