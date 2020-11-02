from drawNA.oxdna import Strand, System
from drawNA.polygons import Edge
import numpy as np
import pandas as pd

edge_1 = Edge(np.array([0.0, 0.0, 0.0]), np.array([5.0, 0.0, 0.0]))
edge_2 = Edge(np.array([5.0, 0.0, 0.0]), np.array([15.0, 0.0, 0.0]))

# Create oxDNA strand
strand_ssDNA = edge_1.strand()
print("ssDNA: ", strand_ssDNA)

strand_dsDNA = edge_2.strand(double=True)
# Create an oxDNA system
oxDNA_system = System(np.array([20.0, 20.0, 20.0]))
strands_to_add = strand_ssDNA + strand_dsDNA
oxDNA_system.add_strands(strands_to_add)
# Print the System class's dataframe
oxDNA_system.write_oxDNA('out')
# Write the system to file
oxDNA_system.dataframe
