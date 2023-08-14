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

use_IR  = False
use_UV  = False
use_IR_spec = False 
use_rho = True  
use_IR_spec_exp = False 


in_rho  = "COMSOL_pressure_scan.txt"
in_IR   = "Ek.dat"
in_IR_spec = "Speck.dat"
in_UV   = "0.1bar_Subt2__0__17-04-30-844.txt"
in_IR_spec_exp = "1.2bar_2.21e-1mbar_1308032U1.txt"

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
IR SPEC EXP FILES [SPECTRUM AS MEASURED BY SPECTROMETER]

Should contain (at least) five columns:
    - wavelength data in nm (col 0)
    - spectral intensity data in arb. units (col 4)

Delimiter: " "

8 lines of header; no comments

"""
if use_IR_spec_exp:

    out_IR_spec_exp = "IRspec_exp.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR_spec_exp), usecols=(0,4), skiprows=8, delimiter=";")
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_IR_spec_exp),arr)

""" 
IR SPEC FILES

Should contain (at least) two columns:
    - wavelength data in nm 
    - spectral intensity data in arb. units 

Delimiter: " "

No headers or comments

"""
if use_IR_spec:

    out_IR_spec = "IRspec.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR_spec), usecols=(0,1))
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_IR_spec),arr)    

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
    #arr[:,0] = scipy.constants.c / (arr[:,0]*1e-9) * 2 * np.pi
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_UV),arr)
""" 
DENSITY FILES

Should contain (at least) two columns:
    - cell spatial data in m; symmetric about zero! 
    - number density in 1/m^3 

Delimiter: " "

Header comments beginning with %

"""
if use_rho: 

    arr = np.loadtxt(os.path.join(in_dir, in_rho), comments="%", delimiter=",")
    arr[:,0] += np.max(arr[:,0])
    arr[:,0] *= 1e-3 

    for i in np.arange(arr.shape[1]-1)+1:
        dens = arr[:,i]
        out_rho="dens_{0:.1f}bar.dat".format(0.1+(i-1)*0.1)
        out_arr =np.empty(shape=(len(dens),2))
        out_arr[:,0]= arr[:,0]
        out_arr[:,1]= dens 
    
        np.savetxt(os.path.join(out_dir, out_rho),out_arr)

