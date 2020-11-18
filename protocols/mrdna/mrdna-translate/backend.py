from imports import *

from softnanotools.logger import Logger
logger = Logger('Manager')

import subprocess
import shutil

from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from drawNA.readers import OXDNAReader
from drawNA.oxdna.strand import generate_helix
from drawNA.oxdna import Nucleotide, Strand, System

class Manager:
    """Container for all functions"""
    def __init__(self, number, spacing, tacoxDNA):
        logger.info('Initialising Manager...')
        self.N = number
        self.spacing = spacing

        # check that path is correct and register as attribute
        self.register_tacoxDNA(tacoxDNA)

        logger.info('Finished setting up Manager')

    def register_tacoxDNA(self, tacoxDNA):
        self._tacoxDNA = tacoxDNA
        logger.info(f'Checking that tacoxDNA exists at {self._tacoxDNA}')
        assert Path(self._tacoxDNA).exists()
        return 

    def tacoxDNA(self, function: str) -> str:
        return str(Path(self._tacoxDNA) / 'src' / f'{function}.py')

    @property
    def box(self):
        """Returns the box size based on the number of replicas 
        and their spacing"""
        _box = np.array([5., 5., 5.]) * self.spacing * self.N
        return _box

    def generate_system(self):
        """Generates an oxDNA system"""
        system = System(self.box)
        return

    def export_PDB(self, system, name):
        """Exports drawNA system to PDB"""
        system.write(name)
        subprocess.check_output([
            self.tacoxDNA("oxDNA_PDB"),
            f'oxdna.{name}.conf',
            f'oxdna.{name}.top',
        ])
        return

    def simulate(self):
        """Executes an mrDNA simulation"""
        return

    def import_PDB(self, name):
        """Imports PDB as drawNA system"""
        subprocess.check_output([
            self.tacoxDNA('PDB_oxDNA'),
            name
        ])
        conf = f'{name}.dat'
        top = f'{name}.top'
        system = OXDNAReader([conf, top]).system
        return

    def process_system(self):
        """Creates clones of systems, translates them and writes to file"""
        return