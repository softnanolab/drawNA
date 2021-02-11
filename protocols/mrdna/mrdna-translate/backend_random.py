import json
import re
import os
import shutil
import subprocess
import glob

from pathlib import Path

import pandas as pd
import numpy as np

from imports import *
import generate_random

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
        overall_length: int, #n_strands: int
        percent_stapled: int, #length: int
        length_stapled: int, #stapled: int
        sim_index: int,
        tacoxDNA: str,
        root: str = 'data',
        name: str = 'main',
        energy: str = 'oxdna.energy',
        oxDNA: str = 'oxDNA'
       # binary: str# = '20_string.txt'
    ):
        self.box = box
        self.number = number
        # self.n_strands = n_strands #new
        # NEW VARIABLES(ALFONSO)
        self.overall_length = int(overall_length)
        self.percent_stapled = percent_stapled
        self.length_stapled = int(length_stapled)

        self.simulation_index = sim_index
        #self.binary_file = binary

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

        self.energy = energy




    def generate_system(self):
        """Generate a simple oxDNA system"""
        system = generate_random.generate_system_all_random(
            [self.box] * 3,
            n_strands=self.n_strands,  #or number
            length=self.length,
            stapled=self.length_stapled,
            percentage=self.percent_stapled
        ) 
        self.system = system
        return system

    def generate_exact_system(self):
        """Generate a simple oxDNA system"""
        system = generate_random.generate_exact_system(
            [self.box] * 3,
            n_strands=self.n_strands,  #or number
            length=self.length,
            stapled=self.length_stapled,
            percentage=self.percent_stapled,
            binary_file='20_string.txt'#self.binary_file
        ) 
        self.system = system
        return system

    def generate_loop_system(self):
        """Generate a simple oxDNA system"""
        system = generate_random.generate_loop_random(
            [self.box] * 3,
            n_strands=self.n_strands,  #or number
            length=self.length,
            stapled=self.length_stapled,
            percentage=self.percent_stapled,
        ) 
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
            CUDA_device = 0,
            backend_precision='mixed',
            steps=10000, #1000000,
            newtonian_steps=103,
            diff_coeff=2.50,
            thermostat='john',
            T=0.1, #0.110 #0.1 represents 300K. However, at room temp equilibration takes very long
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
            print_conf_interval=10000, #10000,
            print_energy_every=1000, #1000,
            trajectory_file=f'{self.file_path}/oxdna.{self.name}.traj',
            lastconf_file=self.configuration
        )

        if stage == 'equilibration':
            pass

        elif stage == 'replication':
            input_file['print_conf_interval'] = 10000#100000
            input_file['steps'] = 10000 * self.number #100000 * self.number

      #  elif stage == 'box_reduction':
      #      input_file['print_conf_interval'] = 100000
      #      input_file['steps'] = 100000 * self.number

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
        if abs(gradient) < 1e-1: #1e-3:
            return True
        else:
            return False

    def reduce_box_size(self):
        """reduction of the box size"""
        step_size = 0.0

        with open(self.configuration, 'w') as conf:
            data = file.readlines()
            line = data[1]
            line = line.replace('b', '').replace('=', '')
            numbers = np.fromstring(line)

            current_box_size = numbers[0]
            new_box_size = current_box_size*(1-step_size)
            data[1] = 'b = {} {} {}'.format(new_box_size, numbers[1], numbers[2])
        return

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
                    self._oxDNA, 
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

    def run_box_reduction(self):
        """reduces the box size step by step to increase pressure on the system"""
        counter = 0 #this counter is only for testing right now. in the future a final pressure should be used for the while loop
        logger.info('Starting box reduction...')
        try:

            subprocess.check_output([
                self._oxDNA, 
                f'{self.root}/oxdna.box_reduction.input'
            ])
        
            while counter < 1:
                counter += 1
                self.reduce_box_size()
                logger.info('Final pressure not reached. Box will be further reduced...')
                subprocess.check_output([
                    self._oxDNA, 
                    f'{self.root}/oxdna.box_reduction.input'
                ])

        except KeyboardInterrupt:
            logger.debug('Caught KeyboardInterrupt successfully')

        return

    def split_trajectory(self):
        """Splits the trajectory into sub-configurations"""
        n_lines = len(self.system.nucleotides) + 3
        # get number of lines in trajectory
        with open(self.file_path + '/' + 'oxdna.main.traj', 'r') as f:
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
                #box/4.0 * ((i // m) % m),
                #box/4.0 * (i % m),
                box/4.0 * (i % m),
                box/4.0 * ((i // m) % m),
                box/4.0 * (i // (m * m)),
            ])

            clone.translate(translation)
            mother.add_strands(clone.strands)

        translation = np.array([
            -(box/2.7),
            -(box/2.7),
            -(box/2.7),
        ])
        mother.translate(translation)
        mother.write_oxDNA('mother')



      #  # try reduction of box size instead of applying pressure
      #  logger.info(f'>>> Running box reduction on mother system')
      #  self.topology = f'oxdna.mother.top'
      #  self.configuration = f'oxdna.mother.conf'

       # self.generate_input_file('box_reduction')

       # self.run_box_reduction()

       # logger.info(f'>>> Box reduction finished!')


        logger.info(f'>>> Exporting mother to PDB')
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


        #SIMULATE
        try: 
        # in case the simulation is failing - skip it so that the other simulations in line can run at least. 
        # read_model from mrdna sometimes gives "warning: reconnecting" errors for no apparent reason. if this happens 
        # also skip to the next simulation.

            #READ PDB
            model = read_model('oxdna.mother.conf.pdb')
            simulate(
                model, 
                'out', 
                coarse_steps=1000000000, 
                fine_steps=100,
                coarse_output_period=1000000,
                fine_output_period=100,#,
            #    temperature=198,
                directory='.'
            )
            shutil.move('out-3.pdb', f'{self.root}/mother.pdb')

        except Exception as e:
            print('This failed, skipping.')
            print('The error was: ')
            print(e)
            print('------')
        return

    def cleanup_files(self):
        """Delete all output files"""
        shutil.move(
            f'oxdna.{self.name}.traj',
            f'{self.root}/oxdna.{self.name}.traj'
        )
        shutil.move(
            f'oxdna.energy',
            f'{self.root}/oxdna.energy'
        )
        for pattern in [
            'out-*',
            'rb_*',
            'last*',
            'energy*'
        ]:

            for p in Path(".").glob(pattern):
                p.unlink()

    def define_multiple_systems(self):
        """
        Create multiple folders for each different simulation system 
        and convert terminal input to generate_system input
        """
        print('In define_multiple_systems')

        bases_needed = self.overall_length * (1/(100.0/self.percent_stapled))
        self.n_strands = int(bases_needed/self.length_stapled)



      #  if self.n_strands == 1:
      #      logger.error(
      #          'n_strands is 1 (number of turns is zero)! it has to be at least 2 in order for the code to work.'
      #          'Decrease -ds or increase -l'
      #          )



        #Create simulation output folder
        final_percentage = round(self.n_strands*self.length_stapled)/self.overall_length
        self.length = int(self.overall_length/self.n_strands)

        #adjust lengths between turns, overall_length and percentage if it doesnt end "6" (physical restriction) 
        length_div = (self.length - 6.0) // 10.45
        self.length = int((length_div * 10.45) + 6.0)
        self.overall_length = self.length * self.n_strands
        final_percentage = round(self.n_strands*self.length_stapled)/self.overall_length
            
        self.file_path = "sim{}_l{}ds{}p{}".format(self.simulation_index,
                                                   self.overall_length,
                                                   self.length_stapled,
                                                   int(final_percentage*100) )# % sim
        # self.sim_path = self.file_path + "/sim"
        # parameters['fpath'] = file_path

        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            print('{} is being deleted'.format(self.file_path))
            shutil.rmtree(self.file_path)
            os.mkdir(self.file_path)
        if self.root:
            self.root = self.file_path + '/' + self.root
            Path(self.root).mkdir(exist_ok=True, parents=True)
        
        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')
        print('\tdone with define_multiple_systems')

    def define_random_system(self):
        """
        Create multiple folders for 5 different simulation systems
        and convert terminal input to generate_system_all_random input
        -> this function is used to randomly distribute blocks of ds staples along the total length of strand
        """

        print('In define_random_systems')

        self.n_strands = 20
        self.length = 400

        if self.percent_stapled != 0:
            final_percentage = 1/(100.0/self.percent_stapled)
        else:
            final_percentage = 0

        #Create simulation output folder
        self.file_path = "sim{}_l{}ds{}p{}".format(self.simulation_index,
                                                   self.overall_length,
                                                   self.length_stapled,
                                                   int(final_percentage*100) )# % sim

        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            print('{} is being deleted'.format(self.file_path))
            shutil.rmtree(self.file_path)
            os.mkdir(self.file_path)
        if self.root:
            self.root = self.file_path + '/' + self.root
            Path(self.root).mkdir(exist_ok=True, parents=True)
        
        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')
        print('\tdone with define_random_system')

    def define_loop_system(self):
        """
        Create multiple folders for 5 different simulation systems
        and convert terminal input to generate_system_all_random input
        -> this function is used to randomly distribute blocks of ds staples along the total length of strand
        """

        print('In define_loop_systems')

        self.n_strands = 20
        self.length = 400 #100 #124

        if self.percent_stapled != 0:
            final_percentage = 1/(100.0/self.percent_stapled)
        else:
            final_percentage = 0

        #Create simulation output folder
        self.file_path = "sim{}_l{}ds{}p{}".format(self.simulation_index,
                                                   self.overall_length,
                                                   self.length_stapled,
                                                   int(final_percentage*100) )# % sim

        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            print('{} is being deleted'.format(self.file_path))
            shutil.rmtree(self.file_path)
            os.mkdir(self.file_path)
        if self.root:
            self.root = self.file_path + '/' + self.root
            Path(self.root).mkdir(exist_ok=True, parents=True)
        
        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')
        print('\tdone with define_random_system')


    def define_random_system_single(self):
        """
        Create multiple folders for 5 different simulation systems
        and convert terminal input to generate_system_all_random input
        -> this function is used to randomly distribute blocks of ds staples along the total length of strand
        THE STRAND DOES NOT HAVE TURNS (single line)
        """

        print('In define_random_systems')

        self.n_strands = 1 
        self.length = 8046 

        if self.percent_stapled != 0:
            final_percentage = 1/(100.0/self.percent_stapled)
        else:
            final_percentage = 0

        #Create simulation output folder 
        self.file_path = "sim{}_l{}ds{}p{}".format(self.simulation_index,
                                                   self.overall_length,
                                                   self.length_stapled,
                                                   int(final_percentage*100))

        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            print('{} is being deleted'.format(self.file_path))
            shutil.rmtree(self.file_path)
            os.mkdir(self.file_path)
        if self.root:
            self.root = self.file_path + '/' + self.root
            Path(self.root).mkdir(exist_ok=True, parents=True)
        
        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')
        print('\tdone with define_random_system')


    def define_exact_system(self, binary_file):
        """
        create a system from a binary file containing exact position of ssDNA and dsDNA
        """

        print('In define_exact_system')

        with open(binary_file, 'r') as binary_system:
            DNA_system = binary_system.read()

        self.n_strands = 20
        self.length = 400
        self.percent_stapled = DNA_system.count('1')/8064

        if self.percent_stapled != 0:
            final_percentage = self.percent_stapled
        else:
            final_percentage = 0
         
        #Create simulation output folder   
        self.file_path = "sim{}_l{}ds{}p{}".format(self.simulation_index,
                                                   self.overall_length,
                                                   self.length_stapled,
                                                   int(final_percentage*100) )
        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            print('{} is being deleted'.format(self.file_path))
            shutil.rmtree(self.file_path)
            os.mkdir(self.file_path)
        if self.root:
            self.root = self.file_path + '/' + self.root
            Path(self.root).mkdir(exist_ok=True, parents=True)
        
        self.topology = f'{self.root}/oxdna.{self.name}.top'
        self.configuration = f'{self.root}/oxdna.{self.name}.conf'
        self.configuration_template = f'{self.root}/oxdna.{self.name}.{{}}.conf'

        logger.debug('Initialised the following filenames:')
        logger.debug(f'\tName         : {self.name}')
        logger.debug(f'\tRoot         : {self.root}')
        logger.debug(f'\tEnergy       : {self.energy}')
        logger.debug(f'\tTopology     : {self.topology}')
        logger.debug(f'\tConfiguration: {self.configuration}')
        logger.debug(f'\tTemplate     : {self.configuration_template.format("X")}')
        print('\tdone with define_exact_system')


    def move_sim_files(self):
        # make sure that the folder exists
        try:
            os.mkdir(self.file_path)
        except FileExistsError:
            pass
        
        temp_files_in_path = glob.glob('./*')
        files_in_path = list()
        for f in temp_files_in_path:
            if 'sim' in f:
                continue
            files_in_path.append(f)
        
            shutil.move(os.path.join('./', f), self.file_path)

#    def save_output_info(self):
#        """
#        Calls vmd - saves particle information (.pdb) from all frames 
#        and generates a movie of the simulation
#        """ 
#        import load-mrdna_copy.tcl as run_vmd
    