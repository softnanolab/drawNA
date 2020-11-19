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

        # set attributes
        self.N = number
        self.spacing = spacing
        self.name = 'main'
        self.root = 'data'
        Path(self.root).mkdir(exists_ok=True)

        # check that path is correct and register as attribute
        self.register_tacoxDNA(tacoxDNA)

        # create a set of empty systems to store the final configurations
        self.mother_system = System(self.box)
        self.child_systems = [System(self.box) for i in range(self.N)]

        # create cache to store temporary objects
        self.cache = {}

        logger.info('Finished setting up Manager')

    def register_tacoxDNA(self, tacoxDNA):
        """Verifies the existence of tacoxDNA and registers it"""
        self._tacoxDNA = tacoxDNA
        logger.info(f'Checking that tacoxDNA exists at {self._tacoxDNA}')
        try:
            assert Path(self._tacoxDNA).exists()
        except AssertionError:
            logger.error(f'tacoxDNA cannot be found at {self._tacoxDNA}!')
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
        """Generate a simple oxDNA system"""
        system = System(self.box)
        strand = generate_helix(40, double=False)[0]
        new = []
        for nt in strand.nucleotides[10:30]:
            new.append(nt.make_across())
        new = Strand(nucleotides=new[::-1])
        system.add_strands([strand, new])    
        return system

    def export_PDB(self, system, name):
        """Exports drawNA system to PDB"""
        logger.info(f'>>> Exporting PDB to {name}')
        system.write_oxDNA(name)
        subprocess.check_output([
            'python',
            self.tacoxDNA("oxDNA_PDB"),
            f'oxdna.{name}.top',
            f'oxdna.{name}.conf',
            '35',
        ])
        logger.info(f'<<< DONE!')
        return f'oxdna.{name}.conf.pdb'

    def oxDNA_simulate(self, **kwargs):
        args = [
            self.oxDNA,
            self.input_file,
        ]
        args += [f'{key}={value}' for key, value in kwargs.items()]
        subprocess.check_output(args)
        return

    def mrdna_simulate(self, pdb_file):
        """Executes an mrDNA simulation"""
        # include imports here to avoid caching errors
        from mrdna.simulate import multiresolution_simulation as simulate
        from mrdna.readers import read_atomic_pdb as read_model    

        # read PDB
        model = read_model(pdb_file)

        # simulate
        simulate(
            model, 
            'out', 
            coarse_steps=1, 
            fine_steps=1,
            coarse_output_period=1,
            fine_output_period=1,
            directory='.')
        self.cache['pdb'] = 'out-3.pdb'
        return

    def import_PDB(self, name):
        """Imports PDB as drawNA system"""
        logger.info(f'>>> Importing PDB from {name}')
        subprocess.check_output([
            'python',
            self.tacoxDNA('PDB_oxDNA'),
            name,
            '35'
        ])
        conf = f'{name}.oxdna'
        top = f'{name}.top'
        system = OXDNAReader([conf, top]).system
        logger.info(f'<<< DONE!')
        return system

    def recenter_system(self, system):
        return system

    def process_system(self):
        """Creates clones of systems, translates them and writes to file"""
        return

    def write_all(self):
        """"""
        return
