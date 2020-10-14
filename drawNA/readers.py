import re

from typing import List

import numpy as np
import pandas as pd
import scipy.spatial.transform

from .oxdna import Nucleotide, Strand, System
from .oxdna.utils import quat_to_exyz

def quat_to_matrix(row: pd.Series) -> np.ndarray:
    """
    Converts a pd.Series which has entries:
        [['qw', 'qi', 'qj', 'qk']]

    from a quaternion to a matrix rotation and returns
    it as a np.ndarray.
    """

    quaternion = np.array([
        row['qw'],
        row['qi'],
        row['qj'],
        row['qk']
    ])

    rotation = scipy.spatial.transform.Rotation(quaternion)

    matrix = rotation.as_matrix()

    return matrix

class Reader:

    columns = [
        'x'   , 'y'   , 'z'  ,
        'a1x' , 'a1y' , 'a1z',
        'a3x' , 'a3y' , 'a3z',
        'vx'  , 'vy'  , 'vz' ,
        'Lx'  , 'Ly'  , 'Lz' ,
        'base', 'before', 'after', 'strand',
    ]

    """
    Low-level class that takes a table of nucleotides and returns
    a `drawNA.oxdna.System` instance. Designed to be subclassed,
    a derived class must pass a dataframe (that can be checked using
    static methods) of nucleotides when initialising, preferably using
    `super()` methods.
    """
    def __init__(
        self, 
        nucleotides: pd.DataFrame, 
        metadata: dict,
    ):

        Reader._check_dataframe(nucleotides)
        Reader._check_metadata(metadata)

        self._nucleotides = nucleotides
        self._metadata = metadata
        
        self._validate_strands()

    @staticmethod
    def _check_metadata(metadata: dict):
        if metadata == None:
            return
        assert True

    @staticmethod
    def _check_dataframe(dataframe: pd.DataFrame):
        assert set(dataframe.columns) == set(Reader.columns)

    def _validate_strands(self):
        assert True

    def _nucleotide_from_series(self, series: pd.Series) -> Nucleotide:
        nuc = Nucleotide(
            series['base'],
            np.array([series['x'], series['y'], series['z']]),
            np.array([series['a1x'], series['a1y'], series['a1z']]),
            np.array([series['a3x'], series['a3y'], series['a3z']]),
            v=np.array([series['vx'], series['vy'], series['vz']]),
            L=np.array([series['Lx'], series['Ly'], series['Lz']]),
        )
        nuc._before = series['before']
        nuc._after = series['after']
        return nuc

    @property
    def strands(self) -> List[Strand]:
        strand_list = []
        strand_counter = True
        i = 0

        # strand_counter will be set to False
        # when an empty strand is found
        while strand_counter:
            i += 1
            temp = self._nucleotides[self._nucleotides['strand']==i]
            strand_counter = temp.any().any()

            # Force break if strand_counter==False
            if not strand_counter:
                break
            nucleotide_list = []
            for index in temp.index:
                nucleotide_list.append(
                    self._nucleotide_from_series(temp.loc[index])
                )
            # TODO: sort nucleotide_list to ensure 
            # before and after are correct
            strand_list.append(Strand(nucleotides=nucleotide_list))
        return strand_list
    
    @property
    def system(self) -> System:
        _system = System(
            self._metadata['box'],
            time=self._metadata.get('time', 0),
            E_pot=self._metadata.get('PE', 0.0),
            E_kin=self._metadata.get('KE', 0.0)
        )
        # reverse list is necessary but cannot remember why
        _system.add_strands(self.strands[::-1])
        return _system

class LAMMPSDataReader(Reader):
    def __init__(self, fname: str):
        self._fname = fname
        super().__init__(self.dataframe, metadata=self.metadata)

    @staticmethod
    def detect_filetypes(fname: str):
        raise NotImplementedError

    def convert_quaternions(self, dataframe: pd.DataFrame) -> pd.DataFrame:
        dataframe['matrix'] = dataframe.apply(
            lambda x: quat_to_matrix(x), 
            axis=1
        )
        for i in range(3):
            for j, comp in enumerate(['x', 'y', 'z']):
                dataframe[f'a{i+1}{comp}'] = \
                    dataframe['matrix'].apply(lambda x: x[i][j])
        

    def read_data(self) -> pd.DataFrame:
        meta = self.metadata
        skip_atoms = meta['skip_atoms']
        skip_velocities = meta['skip_velocities']
        skip_ellipsoids = meta['skip_ellipsoids']
        skip_bonds = meta['skip_bonds']
        n_atoms = meta['n_atoms']
        n_bonds = meta['n_bonds']
        # parse atoms
        atoms = pd.read_csv(
            self._fname, 
            skiprows=skip_atoms,
            delim_whitespace=True,
            header=None,
            nrows=n_atoms,
        ).rename(columns={
            0: 'id',
            1: 'n_base',
            2: 'x',
            3: 'y',
            4: 'z',
            5: 'strand',
            6: 'flag',
            7: 'density'
        })
        atoms['base'] = atoms['n_base'].apply(lambda x: {
                1:'A',
                2:'C',
                3:'G',
                4:'T'
            }[x]
        )
        velocities = pd.read_csv(
            self._fname, 
            skiprows=skip_velocities,
            delim_whitespace=True,
            header=None,
            nrows = n_atoms,
        ).rename(columns={
            0: 'id',
            1: 'vx',
            2: 'vy',
            3: 'vz',
            4: 'Lx',
            5: 'Ly',
            6: 'Lz',
        })
        result = pd.merge(atoms, velocities, on='id')
        ellipsoids = pd.read_csv(
            self._fname, 
            skiprows=skip_ellipsoids,
            delim_whitespace=True,
            header=None,
            nrows = n_atoms,
        ).rename(columns={
            0: 'id',
            1: 'Ix',
            2: 'Iy',
            3: 'Iz',
            4: 'qw',
            5: 'qi',
            6: 'qj',
            7: 'qk',
        })
        self.convert_quaternions(ellipsoids)

        result = pd.merge(result, ellipsoids, on='id')
        bonds = pd.read_csv(
            self._fname, 
            skiprows=skip_bonds,
            delim_whitespace=True,
            header=None,
            nrows = n_bonds,
        ).rename(columns={
            0: 'id',
            1: 'type',
            2: 'atom1',
            3: 'atom2',
        })
        result = result.set_index('id')
        for index in bonds.index:
            i = bonds['atom1'][index]
            j = bonds['atom2'][index]
            result.at[i, 'after'] = int(j)
            result.at[j, 'before'] = int(i)

        result['before'] = result['before'].apply(lambda x: -1 if np.isnan(x) else x)
        result['after'] = result['after'].apply(lambda x: -1 if np.isnan(x) else x)

        return result

    @property
    def dataframe(self) -> pd.DataFrame:
        return self.read_data()[Reader.columns]

    @property
    def metadata(self) -> dict:
        box = np.array([0., 0., 0.])
        with open(self._fname, 'r') as f:
            for i, line in enumerate(f.readlines()):
                if re.findall('Atoms', line):
                    skip_atoms = i + 2
                elif re.findall('Velocities', line):
                    skip_velocities = i + 2
                elif re.findall('Ellipsoids', line):
                    skip_ellipsoids = i + 2
                elif re.findall('Bonds', line):
                    skip_bonds = i + 2
                elif re.findall('xlo xhi', line):
                    box[0] = float(line.split()[1]) - float(line.split()[0])
                elif re.findall('ylo yhi', line):
                    box[0] = float(line.split()[1]) - float(line.split()[0])
                elif re.findall('zlo zhi', line):
                    box[0] = float(line.split()[1]) - float(line.split()[0])
                elif re.findall(r"[0-9]+ atoms", line):
                    n_atoms = float(line.split()[0])
                elif re.findall(r"[0-9]+ bonds", line):
                    n_bonds = float(line.split()[0])

        result = {
            'skip_atoms': skip_atoms,
            'skip_velocities': skip_velocities,
            'skip_ellipsoids': skip_ellipsoids,
            'skip_bonds': skip_bonds,
            'n_atoms': n_atoms,
            'n_bonds': n_bonds,
            'box': box,
        }

        return result

class LAMMPSDumpReader(LAMMPSDataReader):
    """
    LAMMPS DumpFile reader designed specifically for LAMMPS simulations
    that have used the oxDNA or oxDNA2 forcefields, and have data/dump
    files configured as per the literature (i.e. not customised)

    LAMMPSDumpReader subclasses LAMMPSDataReader and is also a subclass
    of Reader, which LAMMPSDataReader also subclasses. This architecture
    is slightly confusing but results in very concise code for this class

    The reason for subclassing LAMMPSDataReader is because the configuration
    data file is needed to generate the bonds, which cannot be stored in a 
    dump file.
    """
    def __init__(self, fnames: List[str]):
        data, self._dump = LAMMPSDumpReader.detect_filetypes(fnames)
        LAMMPSDataReader.__init__(self, data)
        
        self._time = 0
        Reader.__init__(self, self.dataframe, metadata=self.metadata)

    @staticmethod
    def detect_filetypes(fnames: List[str]):
        data = fnames[0]
        dump = fnames[1]
        return data, dump

    @property
    def dump(self) -> str:
        return self._dump

    def update_dump(self, fname: str):
        self._dump = fname
        Reader.__init__(self, self.dataframe, metadata=self.metadata)

    def read_dump(self):
        fname = self.dump
        with open(fname, 'r') as f:
            lines = f.readlines()
        for i, line in enumerate(lines):
            if re.findall(r"TIMESTEP", line):
                self._time = int(lines[i+1])
            if re.findall(r"ATOMS", line):
                skip_atoms = i+1
                columns = line.split()[2:]

        dump = pd.read_csv(
            fname, 
            skiprows=skip_atoms, 
            nrows=self.metadata['n_atoms'],
            header=None,
            names=columns,
            delim_whitespace=True,
        ).rename(columns={
            'c_quat[1]': 'qw',
            'c_quat[2]': 'qi',
            'c_quat[3]': 'qj',
            'c_quat[4]': 'qk',
            'tqx':'Lx',
            'tqy':'Ly',
            'tqz':'Lz',
        })

        self.convert_quaternions(dump) 

        return dump

    @property
    def dataframe(self):
        result = pd.DataFrame()
        data = self.read_data()
        dump = self.read_dump()
        dump = pd.merge(
            dump,
            data.reset_index().rename(
                columns={'index':'id'}
            )[['id', 'before', 'after', 'strand', 'base']],
            on='id'
        )
        result = dump[Reader.columns]
        return result

class OXDNAReader(Reader):
    def __init__(self, fnames: List[str]):
        conf, top = OXDNAReader.detect_filetypes(fnames)
        self._metadata = {}
        self._topology = self.read_topology(top)
        self._configuration = self.read_configuration(conf)

        super().__init__(self.dataframe, metadata=self.metadata)

    @staticmethod
    def detect_filetypes(fnames: List[str]):
        conf = fnames[0]
        top = fnames[1]
        return conf, top

    def read_configuration(self, fname: str) -> pd.DataFrame:
        with open(fname, 'r') as f:
            self._metadata['time'] = f.readline().split('=')[-1].strip()
            self._metadata['box'] = np.array(f.readline().split('=')[-1].strip().split())
            _energies = f.readline().split('=')[-1].strip().split()
            self._metadata['PE'] = float(_energies[0])
            self._metadata['KE'] = float(_energies[1])
            data = pd.read_csv(
                f,
                delim_whitespace=True,
                header=None,
                nrows=self._metadata['n_nucleotides'],
            )
            data = data.rename(columns=dict(zip(range(len(Reader.columns)), Reader.columns)))
            
        return data

    def read_topology(self, fname: str) -> pd.DataFrame:
        with open(fname, 'r') as f:
            n_nucleotides, n_strands = [int(i) for i in f.readline().split()]
            self._metadata['n_nucleotides'] = n_nucleotides
            self._metadata['n_strands'] = n_strands
            data = pd.read_csv(fname, delim_whitespace=True, header=None, skiprows=1)
            data = data.rename(columns={
                0: 'strand',
                1: 'base',
                2: 'before',
                3: 'after',
            })
        return data

    @property
    def dataframe(self) -> pd.DataFrame:
        result = pd.concat([self._configuration, self._topology], axis=1)
        return result

    @property
    def metadata(self) -> dict:
        return self._metadata
