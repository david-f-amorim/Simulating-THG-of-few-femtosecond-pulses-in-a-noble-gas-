import numpy as np
import os
import scipy.constants 

# # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#
#          EK READER
#          -------------
#               reads in a "Ek.dat" file (from FROG) and 
#               converts to "pulse.dat" file. 
#
#          Expected input file:
#               - filename: "Ek.dat"
#               - file directory: "./input"
#               - delimiter: " "
#               - data type: float (with decimal points, not commas!)
#               - no headers, comments, etc. 
#               - first col: t (fs)
#               - second col: temporal intensity (arb. units)
#
#         Generated output file:
#               - filename: "pulse.dat"
#               - file directory: "./input"
#               - delimiter: " "
#               - data type: float (with decimal points, not commas!)
#               - no headers, comments, etc. 
#               - first col: t (fs)
#               - second col: temporal intensity (arb. units)
#
# # # # # # # # # # # # # # # # # # # # # # # # # # # # #

dir      = "input"
in_name  = "Ek.dat"
out_name = "pulse.dat"

arr = np.loadtxt(os.path.join(dir, in_name), usecols=(0,1))
arr[:,0] *= 1e-15

np.savetxt(os.path.join(dir, out_name),arr)