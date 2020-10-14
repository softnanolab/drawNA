#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 18:25:44 2020

@author: louiserosset
"""

from drawNA.oxdna import System, Strand, Nucleotide
from drawNA.oxdna.system import oxDNA_string
import os
import pathlib
import numpy as np

class SystemWithFolder(System):
    """
    oxDNA simulation system, a container for strands which
    are composed of nucleotides.

    Parameters:
        box - the box size of the system
        time - Time of the system
        E_pot - Potential energy
        E_kin - Kinetic energy
    """
    
    def __init__(self, box: np.ndarray, time: int = 0, E_pot: float = 0.0, E_kin: float = 0.0):
        super().__init__(box, time, E_pot, E_kin)


    def write_oxDNA_folder(self, sim_nb, tot_l, ds_l, ss_l, prefix: str = "out"):
        """
        Creates one folder to contain the configuration file and topology 
        file required to run a simulation using oxDNA and also outputs an 
        info file with the total length, ds_portion and ss_portion lengths.
        Parameters:
            sim_nb : integer
            prefix ('out') : prefix to output files
        """
        pathname = "oxDNA_sims/sim" + "{:.0f}".format(sim_nb)
        filename1 = f"oxdna.{prefix}.conf"
        filename2 = f"oxnda.{prefix}.top"
        
        pathlib.Path(pathname).mkdir(parents=True, exist_ok=True) 
        
        with open(os.path.join(pathname, filename1), "w") as f:
            f.write(f"t = {self.time}\n")
            f.write(f"b = {self.box[0]} {self.box[1]} {self.box[2]}\n")
            f.write(f"E = {self.E_pot} {self.E_kin} {self.E_tot}\n")
            f.write(oxDNA_string(self.configuration))

        with open(os.path.join(pathname, filename2), "w") as f:
            f.write(f"{len(self.nucleotides)} {len(self.strands)}\n")
            f.write(oxDNA_string(self.topology))

        with open(os.path.join(pathname,"sim_info"), "w") as f:
            f.write("tot length is {:.0f}\n".format(tot_l))
            f.write("ds length is {:.0f}\n".format(ds_l))
            f.write("ss length is {:.0f}\n".format(ss_l))