import numpy as np
import os
import scipy.constants 

"""
 file_prepare.py: 

        Reads in data files containing information about the input (IR) pulse,
        the gas density distribution, and the output (UV) pulse. Writes the
        relevant content to appropriate files, which can then be used in 
        "THG_sim_main.jl".
"""

in_dir  = "raw_input" 
out_dir = "input" 

use_IR  = True
use_UV  = False 
use_rho = False  

in_IR   = "Ek.dat"
in_UV   = "PUT_FILENAME_HERE"
in_rho  = "PUT_FILENAME_HERE"

""" 
IR PULSE FILES

Should contain (at least) two columns:
    - timing data in fs 
    - temporal intensity data in arb. units 

Delimiter: " "

No headers or comments

"""
if use_IR:

    out_IR = "IRpulse.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR), usecols=(0,1))
    arr[:,0] *= 1e-15

    np.savetxt(os.path.join(out_dir, out_IR),arr)

""" 
UV PULSE FILES

Should contain (at least) two columns:
    - wavelengths in nm 
    - spectral intensity data in arb. units 

Delimiter: " "

Contains 14-line header block 

"""
if use_UV: 

    out_UV = "UVpulse.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_UV), usecols=(0,1), skiprows=14)
    arr[:,0] = scipy.constants.c / (arr[:,0]*1e-9) * 2 * np.pi

    np.savetxt(os.path.join(out_dir, out_UV),arr)
""" 
DENSITY FILES

Should contain (at least) two columns:
    - cell spatial data on the interval [0,L] in m 
    - number density in 1/m^3 

Delimiter: " "

No headers or comments

"""
if use_rho: 

    out_rho = "dens.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_rho), usecols=(0,1))

    np.savetxt(os.path.join(out_dir, out_rho), arr)

