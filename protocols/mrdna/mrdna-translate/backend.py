import json
import re

from pathlib import Path

from imports import *

from softnanotools.logger import Logger
logger = Logger('Manager')

from drawNA.oxdna.strand import generate_helix
from drawNA.oxdna import Nucleotide, Strand, System

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
        name: str = 'main'
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
            
        self.tacoxDNA = tacoxDNA

        # initialise filenames
        self.name = name
        self.root = root
        if self.root:
            Path(self.root).mkdir(exist_ok=True, parents=True)


        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')

    def generate_system(self) -> System:
        """Creates an oxDNA system"""
        system = System([self.box] * 3)
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
            r0=self.box,
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
        input_file = dict(
            backend='GPU',
            backend_precision='double',
            steps=10,
            newtonian_steps=103,
            diff_coeff=2.50,
            thermostat='john',
            T=0.110,
            dt=0.005,
            verlet_skin=0.05,
            topology=self.topology,
            configuration=self.configuration,
            refresh_vel=1,
            restart_step_counter=1,
            energy_file='oxdna.energy',
            time_scale='linear',
            external_forces=True,
            external_forces_file=self.forces,
            print_conf_interval=10,
            print_energy_every=10,
        )
        string = oxDNA_string(input_file)
        string = re.sub(r"[\{|\}]", "", string).replace("  ", "")
        with open(f'{self.root}/oxdna.input', 'w') as f:
            f.write(string[1:])
        return

    def run_equilibration(self):
        return

    def run_replication(self):
        return
