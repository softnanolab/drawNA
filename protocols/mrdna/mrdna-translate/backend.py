import json
import re
import os
import shutil
import subprocess

from pathlib import Path

import pandas as pd
import numpy as np

from imports import *
import generate

from softnanotools.logger import Logger
logger = Logger('Manager')

from drawNA.oxdna import Nucleotide, Strand, System
from drawNA.readers import OXDNAReader

def oxDNA_string(dictionary) -> str:
    string = json.dumps(dictionary, indent=2)
    string = string.replace('"', '' 
        ).replace(":", " ="
        ).replace(',\n', '\n'
    )
    return string

class Manager:
    def __init__(
        self, 
        box: float, 
        number: int, 
        spacing: float,
        tacoxDNA: str,
        root: str = 'data',
        name: str = 'main',
        energy: str = 'oxdna.energy',
        oxDNA: str = 'oxDNA'
    ):
        self.box = box
        self.number = number
        self.spacing = spacing
        
        # check tacoxDNA is installed properly
        if tacoxDNA == None:
            logger.error(
                'tacoxDNA has not been set! It can be set by doing either:\n'
                '\t1. export TACOXDNA=/path/to/tacoxDNA\n'
                '\t2. python main.py ... --tacoxDNA /path/to/tacoxDNA'
            )
            
        self._tacoxDNA = tacoxDNA
        self._oxDNA = oxDNA

        # initialise filenames
        self.name = name
        self.root = root
        if self.root:
            Path(self.root).mkdir(exist_ok=True, parents=True)


        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'
        self.energy = energy

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')

    def generate_system(self):
        """Generate a simple oxDNA system"""
        system = generate.generate_system([self.box] * 3) 
        self.system = system
        return system

    def generate_forces(self) -> str:
        """Generates a file that ensures a 
        spherical inclusion potential is applied"""
        forces_filename = f'{self.root}/oxdna.forces'
        forces = dict(
            type='sphere',
            particle=-1,
            center="0., 0., 0.",
            r0=self.box/2,
            stiff=5.0,
        )
        string = oxDNA_string(forces)
        with open(forces_filename, 'w') as f:
            f.write(string)

        logger.debug(f'The following forces have been written to {forces_filename}')
        logger.debug(string)

        self.forces = forces_filename
        return forces_filename

    def generate_input_file(self, stage):
        """Outputs oxDNA input files"""
        # create input file dictionary
        input_file = dict(
            sim_type='MD',
            backend='CUDA',
            backend_precision='mixed',
            steps=100000,
            newtonian_steps=103,
            diff_coeff=2.50,
            thermostat='john',
            T=0.110,
            dt=0.005,
            verlet_skin=0.05,
            topology=self.topology,
            conf_file=self.configuration,
            refresh_vel=1,
            restart_step_counter=1,
            energy_file=self.energy,
            time_scale='linear',
            external_forces=True,
            external_forces_file=self.forces,
            print_conf_interval=100000,
            print_energy_every=1000,
            trajectory_file='oxdna.main.traj',
            lastconf_file=self.configuration
        )

        if stage == 'equilibration':
            pass

        elif stage == 'replication':
            input_file['print_conf_interval'] = 10000
            input_file['steps'] = 10000 * self.number

        # format string
        string = oxDNA_string(input_file)
        string = re.sub(r"[\{|\}]", "", string).replace("  ", "")

        # write to file
        with open(f'{self.root}/oxdna.{stage}.input', 'w') as f:
            f.write(string[1:])

        return

    def check_energy(self) -> bool:
        """Checks the energy of the system to see
        if equilibration has been achieved"""
        data = pd.read_csv(
            self.energy, 
            delim_whitespace=True, 
            header=None
        )
        gradient = np.gradient(data[3][-50:]).mean()
        logger.debug(f'Checking {self.energy}, gradient={gradient}')
        if abs(gradient) < 1e-3:
            return True
        else:
            return False

    def run_equilibration(self):
        """Executes the equilibration stage"""
        self.system.write_oxDNA(self.name)
        shutil.move(
            f'./oxdna.{self.name}.conf',
            self.configuration
        )
        shutil.move(
            f'./oxdna.{self.name}.top',
            self.topology
        )
        
        logger.info('Calling oxDNA...')
        try:
            subprocess.check_output([
                self._oxDNA, 
                f'{self.root}/oxdna.equilibration.input'
            ])
        

            while not self.check_energy():
                logger.info('Checking energy...')
                subprocess.check_output([
                    oxDNA, 
                    f'{self.root}/oxdna.equilibration.input'
                ])
        except KeyboardInterrupt:
            logger.debug('Caught KeyboardInterrupt successfully')

        return

    def run_replication(self):
        """Executes the replication stage"""
        subprocess.check_output([
            self._oxDNA, 
            f'{self.root}/oxdna.replication.input'
        ])
        return

    def split_trajectory(self):
        """Splits the trajectory into sub-configurations"""
        n_lines = len(self.system.nucleotides) + 3
        # get number of lines in trajectory
        with open('oxdna.main.traj', 'r') as f:
            data = f.readlines().copy()

            for i in range(self.number):
                fname = self.configuration_template.format(i)
                logger.debug(f'Writing replica to {fname}')
                string = ''.join(data[n_lines*i : n_lines*(i+1)])
                with open(fname, 'w') as fout:
                    fout.write(string)
        return

    def tacoxDNA(self, function: str) -> str:
        """Returns string of tacoxDNA executables"""
        return str(Path(self._tacoxDNA) / 'src' / f'{function}.py')

    def create_mother_system(self):
        """Creates a large mother system"""
        m = np.ceil(np.cbrt(self.number)) 
        box = m * self.box
        mother = System([box, box, box])
        for i in range(self.number):
            clone = OXDNAReader([
                self.configuration_template.format(i), 
                self.topology
            ]).system

            translation = np.array([
                box/2 * (i % m),
                box/2 * ((i // m) % m),
                box/2 * (i // (m * m)),
            ])

            clone.translate(translation)
            mother.add_strands(clone.strands)
        mother.write_oxDNA('mother')

        logger.info(f'>>> Exporting PDB to mother')
        subprocess.check_output([
            'python',
            self.tacoxDNA("oxDNA_PDB"),
            f'oxdna.mother.top',
            f'oxdna.mother.conf',
            '35',
        ])
        logger.info(f'<<< DONE!')
        return

    def mrdna_simulate(self):
        """Executes an mrDNA simulation"""
        # include imports here to avoid caching errors
        from mrdna.simulate import multiresolution_simulation as simulate
        from mrdna.readers import read_atomic_pdb as read_model    

        # read PDB
        model = read_model('oxdna.mother.conf.pdb')

        # simulate
        simulate(
            model, 
            'out', 
            coarse_steps=100, 
            fine_steps=100,
            coarse_output_period=1,
            fine_output_period=1,
            directory='.'
        )
        shutil.move('out-3.pdb', f'{self.root}/mother.pdb')
        return

    def cleanup_files(self):
        """Delete all output files"""
        shutil.move(
            f'oxdna.{self.name}.traj'
            f'{self.root}/oxdna.{self.name}.traj'
        )
        shutil.move(
            f'oxdna.{self.name}.energy'
            f'{self.root}/oxdna.{self.name}.energy'
        )
        for pattern in [
            'out*',
            'rb_*',
            'last*',
            'energy*'
        ]:

            for p in Path(".").glob(pattern):
                p.unlink()