import uclchem
import numpy as np
import matplotlib.pyplot as plt

################################################
# User Inputs Go Here
################################################

species_name = "#CO"
result_file = "examples/test-output/phase1-full.dat"
output="analysis.dat"

################################################
uclchem.analysis.analysis(species_name, result_file, output)