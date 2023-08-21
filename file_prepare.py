import numpy as np
import os

"""
 file_prepare.py: 

        Reads in raw data files containing information about the input (IR) pulse,
        the gas density distribution, and/or the output (UV) pulse. Writes the
        relevant content to appropriate files, which can then be used in 
        "THG_sim_main.jl".
"""


# ----------------- QUICK SETTINGS -------------------------------

use_IR          = False      # if true, IR time-intensity data file is processed
use_UV          = True       # if true, UV spectrum data file is processed
use_IR_spec     = False      # if true, IR spectrum data file is processed 
use_rho         = False      # if true, density data file is processed 
use_IR_spec_exp = False      # if true, a second IR spectrum data file is processed (NOTE: mainly legacy)

# ----------------- FILE NAMES AND PATHS -------------------------------

in_dir  = "raw_input"    # name of input directory containing the raw data
out_dir = "input"        # name of input directory to which processed data is saved  (NOTE: should be the same as the "THG_sim_main.jl" input directory) 

in_rho         = "COMSOL_pressure_scan.txt"             # name of file containing density data 
in_IR          = "Ek.dat"                               # name of file containing IR time-intensity information
in_IR_spec     = "Speck.dat"                            # name of file containing IR spectrum data 
in_UV          = "1.0bar_Subt2__0__17-08-23-927.txt"    # name of file containing UV spectrum data 
in_IR_spec_exp = "1.2bar_2.21e-1mbar_1308032U1.txt"     # name of file containing alternative IR spectrum data 

# ----------------- FILE PROCESSING -------------------------------

# * * * IR time-intensity files * * * 
#
#   - File is expected to contain (at least) two columns:
#       - first column:  time in fs 
#       - second column: intensity in arb. units 
#
#   - Expected delimiter: " "
#
#   - Should contain no headers or comments
#
#   - Processing removes unnecessary columns and 
#     converts to seconds 

if use_IR:

    out_IR = "IRpulse.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR), usecols=(0,1))
    arr[:,0] *= 1e-15

    np.savetxt(os.path.join(out_dir, out_IR),arr)

# * * * IR spectrum files (as measured by spectrometre) * * * 
#
#   - Should contain (at least) five columns:
#       - first column: wavelength data in nm 
#       - fifth column: intensity data in arb. units 
#   - Expected delimiter: " "
#
#   - Should contain 8 lines of header; no comments
#
#   - Processing removes unnecessary columns and 
#     header and converts to metres 

if use_IR_spec_exp:

    out_IR_spec_exp = "IRspec_exp.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR_spec_exp), usecols=(0,4), skiprows=8, delimiter=";")
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_IR_spec_exp),arr)

# * * * IR spectrum files (as produced by FROG) * * * 
#
#   -Should contain (at least) two columns:
#       - first column: wavelength data in nm 
#       - second column: intensity data in arb. units 
#   - Expected delimiter: " "
#
#   - Should contain no header or comments
#
#   - Processing removes unnecessary columns and 
#     header and converts to metres 

if use_IR_spec:

    out_IR_spec = "IRspec.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_IR_spec), usecols=(0,1))
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_IR_spec),arr)    

# * * * UV spectrum files (as produced by FROG) * * * 
#
#   -Should contain (at least) two columns:
#       - first column: wavelength data in nm 
#       - second column: intensity data in arb. units 
#   - Expected delimiter: " "
#
#   - Should contain 14-line header; no comments
#
#   - Processing removes unnecessary columns and 
#     header and converts to metres 

if use_UV: 

    out_UV = "UVpulse.dat"

    arr = np.loadtxt(os.path.join(in_dir, in_UV), usecols=(0,1), skiprows=14)
    arr[:,0] *= 1e-9

    np.savetxt(os.path.join(out_dir, out_UV),arr)

# * * * density files * * * 
#
#   -Should contain (at least) two columns:
#       - first column: position data in mm, centred on zero 
#       - second column: number density in 1/m^3 at central pressure 0.1 bar 
#       - nth column: number density in 1/m^3 at central pressure (0.1+ 0.1 [n-2])bar      
#   - Expected delimiter: ","
#
#   - Contains comments starting with "%"
#
#   - Processing removes centring on zero, converts to metre, and 
#     saves the density profile of each pressure value in a separate file

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

