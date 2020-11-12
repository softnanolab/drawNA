"""
oxDNA system class - contains tools to manage, read and write
oxDNA simulation configurations.
"""
import re
from typing import List

import numpy as np
import pandas as pd
import pathlib
import os

from drawNA.oxdna import Strand, Nucleotide

CONFIGURATION_COLUMNS = ["position", "a1", "a3", "v", "L"]
TOPOLOGY_COLUMNS = ["strand", "base", "3p", "5p"]

LMP_COL_ATOMS = ["id", "type", "position", "molecule", "flag", "density"]
LMP_COL_VELOCITIES = ["id", "v", "L"]
LMP_COL_ELLIPSOIDS = ["id", "shape", "quaternion"]


def oxDNA_string(dataframe: pd.DataFrame) -> str:
    """
    Formats the dataframes needed for writing the topology
    and configuration dataframes to the appropriate
    file format
    """
    output = dataframe.to_string(
        header=False,
        index=False,
        justify="left",
        # format to ensure 4 decimal places
        formatters={
            "position": lambda x: [f"{i:.4f}" for i in x],
            "a1": lambda x: [f"{i:.4f}" for i in x],
            "a3": lambda x: [f"{i:.4f}" for i in x],
            "v": lambda x: [f"{i:.4f}" for i in x],
            "L": lambda x: [f"{i:.4f}" for i in x],
        },
    )
    # remove all excess symbols and whitespace
    output = re.sub(r"\[|\]|\'|\`|\,", "", output)
    # there are a combination of triple & double spaces
    # hence substitute all multi-spaces with one space
    output = re.sub(r" +", " ", output)
    output = output.strip()
    output = output.replace("\n ", "\n")
    return output


def lammps_string(dataframe: pd.DataFrame) -> str:
    """
    Formats the dataframes needed for writing the
    lammps dataframe to the appropriate file format.
    """
    output = dataframe.to_string(
        header=False,
        index=False,
        justify="left",
        # format to ensure 4 decimal places
        formatters={
            "position": lambda x: [f"{i:.10f}" for i in x],
            "shape": lambda x: [f"{i:.10f}" for i in x],
            "quaternion": lambda x: [f"{i:.10f}" for i in x],
            "v": lambda x: [f"{i:.10f}" for i in x],
            "L": lambda x: [f"{i:.10f}" for i in x],
        },
    )
    # remove all excess symbols and whitespace
    output = re.sub(r"\[|\]|\'|\`|\,", "", output)
    # there are a combination of triple & double spaces
    # hence substitute all multi-spaces with one space
    output = re.sub(r" +", " ", output)
    output = output.strip()
    output = output.replace("\n ", "\n")
    return output


class System:
    """
    oxDNA simulation system, a container for strands which
    are composed of nucleotides.

    Parameters:
        box - the box size of the system
        time - Time of the system
        E_pot - Potential energy
        E_kin - Kinetic energy
    """

    def __init__(
        self, box: np.ndarray, time: int = 0, E_pot: float = 0.0, E_kin: float = 0.0
    ):

        self.box = box
        self.time = float(time)
        self.E_pot = E_pot
        self.E_kin = E_kin

        self._strands = []

    def __repr__(self) -> str:
        return f"oxDNASystem[strands: {len(self._strands)}, nucleotides: {len(self.nucleotides)}]"

    @property
    def E_tot(self):
        return self.E_pot + self.E_kin

    @property
    def bonds(self) -> int:
        """
        Returns the total number of bonds in the system,
        used for writing LAMMPS configuration data files.
        """
        result = pd.concat(i.bonds for i in self.strands)
        result["id"] = [i + 1 for i in range(len(result))]
        return result[["id", "type", "atom_1", "atom_2"]]

    @property
    def lammps(self) -> List[pd.DataFrame]:
        """
        Returns a DataFrame containing the
        information needed to write a LAMMPS configuration
        data file
        """
        result = pd.concat(i.lammps for i in self.strands)
        return result

    @property
    def strands(self) -> list:

        # used to set Strand._nucleotide_shift
        shift = 0
        for i, strand in enumerate(self._strands):
            strand._nucleotide_shift = shift
            # required for oxDNA .top format
            strand.index = i + 1
            shift += len(strand)
        return self._strands

    @property
    def nucleotides(self) -> list:
        result = []
        for strand in self.strands:
            result += strand._nucleotides
        return result

    @property
    def dataframe(self) -> pd.DataFrame:
        return pd.concat([i.dataframe for i in self.strands]).reset_index(drop=True)

    @property
    def configuration(self) -> pd.DataFrame:
        return self.dataframe[CONFIGURATION_COLUMNS]

    @property
    def topology(self) -> pd.DataFrame:
        return self.dataframe[TOPOLOGY_COLUMNS]

    def write_oxDNA(self, prefix: str = "out", root: str = "."):
        """
        Writes two files oxdna.*.conf and oxdna.*.top for the
        configuration file and topology file required
        to run a simulation using oxDNA

        Parameters:
            prefix ('out') : prefix to output files
        """
        with open(f"{root}/oxdna.{prefix}.conf", "w") as f:
            f.write(f"t = {self.time}\n")
            f.write(f"b = {self.box[0]} {self.box[1]} {self.box[2]}\n")
            f.write(f"E = {self.E_pot} {self.E_kin} {self.E_tot}\n")
            f.write(oxDNA_string(self.configuration))

        with open(f"{root}/oxdna.{prefix}.top", "w") as f:
            f.write(f"{len(self.nucleotides)} {len(self.strands)}\n")
            f.write(oxDNA_string(self.topology))

    def write_LAMMPS(self, prefix: str = "out", root: str = "."):
        """
        Writes lammps.*.conf which is a configuration data
        file needed to run a lammps simulation.
        
        Parameters;
            prefix ('out') : prefix to output file
        """

        with open(f"{root}/lammps.{prefix}.conf", "w") as f:
            f.write(f"# LAMMPS data file\n")
            f.write(f"{len(self.nucleotides)} atoms\n")
            f.write(f"{len(self.nucleotides)} ellipsoids\n")
            f.write(f"{len(self.bonds)} bonds\n\n")
            f.write(f"4 atom types\n")
            f.write(f"1 bond types\n\n")
            f.write(f"0.0 {self.box[0]} xlo xhi\n")
            f.write(f"0.0 {self.box[1]} ylo yhi\n")
            f.write(f"0.0 {self.box[2]} zlo zhi\n\n")

            f.write(f"Masses\n\n")
            f.write(f"1 3.1575\n")
            f.write(f"2 3.1575\n")
            f.write(f"3 3.1575\n")
            f.write(f"4 3.1575\n")

            f.write("\nAtoms\n\n")
            f.write(lammps_string(self.lammps[LMP_COL_ATOMS]))
            f.write("\n\nVelocities\n\n")
            f.write(lammps_string(self.lammps[LMP_COL_VELOCITIES]))
            f.write("\n\nEllipsoids\n\n")
            f.write(lammps_string(self.lammps[LMP_COL_ELLIPSOIDS]))
            f.write("\n\nBonds\n\n")
            f.write(lammps_string(self.bonds))
            f.write("\n")

    def add_strand(self, addition: Strand, index: int = None):
        """
        Method to add strand(s) to the system

        Parameters:
            addition - accepted as Strand objects or a List of Strands
            index (default = None) - Strand will append to current system,
            otherwise Strand inserted at location given
                
        """

        try:
            assert isinstance(addition, Strand)
        except TypeError:
            raise TypeError(f"addition must be Strand but is {type(addition)}")

        try:
            if index == None:
                self._strands.append(addition.copy())
            else:
                self._strands.insert(index, addition.copy())
        except TypeError:
            raise TypeError("Index must an an integer")

    def add_strands(self, strand_obj: (list or dict) = None, index: int = None):
        """
        Add multiple strands to the system. Use a list of strands with
        an index indicating where they start, or a dictionary where
        the key is an integer and the value is a strand. 
        
        strand_lists are added in reverse order 
        to preserve original list order

        strand_dicts are added in normal order
        to preserve dictionary keys

        Parameters:
            - strand_obj (None) : list or dict of strands
            - index (None) : index to start adding list
        """
        if isinstance(strand_obj, list):
            for strand in strand_obj[::-1]:
                self.add_strand(strand, index)

        elif isinstance(strand_obj, dict):
            for index, strand in sorted(strand_obj.items()):
                self.add_strand(strand, index)

        else:
            raise TypeError(
                "add_strands() requires ONE of a list or dictionary of strands"
            )

    def transform(self, matrix: np.ndarray):
        for strand in self._strands:
            for nucleotide in strand._nucleotides:
                nucleotide.transform(matrix)
        return

    def translate(self, translation_vector: np.ndarray):
        for strand in self._strands:
            for nucleotide in strand._nucleotides:
                nucleotide.translate(translation_vector)
        return

    def rotate(self, rotator: np.ndarray):
        for strand in self._strands:
            for nucleotide in strand._nucleotides:
                nucleotide.rotate(rotator)
        return
