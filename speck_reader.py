import numpy as np
import os
import scipy.constants 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#          SPECK READER
#          -------------
#               reads in a "Speck.dat" file (from FROG) and 
#               converts to "DataField.dat" file. 
#
#          Expected input file:
#               - filename: "Speck.dat"
#               - file directory: "./input"
#               - delimiter: " "
#               - data type: float (with decimal points, not commas!)
#               - no headers, comments, etc. 
#               - first col: lambda (nm)
#               - second col: spectral intensity (arb. units)
#               - third col: unwrapped spectral phase (radians)
#
#         Generated output file:
#               - filename: "DataField.dat"
#               - file directory: "./input"
#               - delimiter: " "
#               - data type: float (with decimal points, not commas!)
#               - no headers, comments, etc. 
#               - first col: frequency (Hz)
#               - second col: spectral intensity (arb. units)
#               - third col: unwrapped spectral phase (radians)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dir      = "input"
in_name  = "Speck.dat"
out_name = "DataField.dat"

arr = np.loadtxt(os.path.join(dir, in_name), usecols=(0,1,2))
arr[:,0] = scipy.constants.c  / (arr[:,0]* 1e-9)

np.savetxt(os.path.join(dir, out_name),arr)