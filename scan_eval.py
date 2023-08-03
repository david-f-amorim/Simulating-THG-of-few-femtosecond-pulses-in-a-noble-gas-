import numpy as np
import matplotlib.pyplot as plt
import os  

"""
 scan_eval.py: 

        Processing and visualisation of pressure scan output. 
"""


scan_dir = "output\scan_150.0mW_Ar_0.0rad_f_ion" 

energy     = True 
efficiency = True 
spectrum   = True 


# ------------- GET PRESSURES, ENERGIES, EFFICIENCIES -------------
arr    = np.loadtxt(os.path.join(scan_dir,"energy_efficiency.txt"))
arr    = arr[arr[:, 0].argsort()]
p_arr  = arr[:,0] 
en_arr = arr[:,1]
ef_arr = arr[:,2]

# ------------- GET SPECTRA ---------------------------------------
N = len(p_arr)
data = np.empty(shape=(N,3), dtype="object")

for i in np.arange(N):

    tmp_arr = np.loadtxt(os.path.join(scan_dir,str(p_arr[i])+" bar.dat"))
    lam = tmp_arr[:,0]
    I = tmp_arr[:,1]

    data[i,0] = p_arr[i]
    data[i,1] = lam 
    data[i,2] = I 

# ---------- GET PARAMS -----------------------------------------

params = np.loadtxt(os.path.join(scan_dir, "params.txt"),dtype="str", delimiter="=", comments="#")

gas  = params[0,1]
phi = float(params[6,1])
beam_en = float(params[7,1])
kerr = params[9,1]
ion = params[10,1]

# ---------- PLOT ENERGY --------------------------------------
peak_energy = np.max(en_arr)
p_at_peak = p_arr[np.where(en_arr == peak_energy )][0]

if energy:
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title(gas+"; {0}mW; {1:.2f}rad; ".format(beam_en*1e6, phi)+kerr+"; ion="+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, en_arr*1e9, color="blue")
    plt.scatter(p_arr, en_arr*1e9, color="blue", label="Peak: {0:.1f}nJ at {1}bar".format(peak_energy*1e9, p_at_peak))
    plt.legend()

    plt.savefig(os.path.join(scan_dir,"energies.png"),dpi=1000)
    plt.show()

# ---------- PLOT EFFICIENCY --------------------------------------
peak_eff = np.max(ef_arr)
p_at_peak_eff = p_arr[np.where(ef_arr == peak_eff )][0]

if efficiency:
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title(gas+"; {0}mW; {1:.2f}rad; ".format(beam_en*1e6, phi)+kerr+"; ion="+ion, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, ef_arr*1e2, color="blue")
    plt.scatter(p_arr, ef_arr*1e2, color="blue", label="Peak: {0:.2f}% at {1}bar".format(peak_eff*1e2, p_at_peak_eff))
    plt.legend()

    plt.savefig(os.path.join(scan_dir,"efficiencies.png"),dpi=1000)
    plt.show()

# ---------- PLOT SPECTRA --------------------------------------
if spectrum:
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV spectra", fontsize=16)
    plt.title(gas+"; {0}mW; {1:.2f}rad; ".format(beam_en*1e6, phi)+kerr+"; ion="+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Wavelength (nm)")

    cmap = plt.get_cmap("viridis")
    cidx = p_arr / np.max(p_arr) 

    for i in np.arange(N):
        plt.plot(data[i,1]*1e9, data[i,2], color=cmap(cidx[i]), label="{0}bar".format(data[i,0]))
    plt.legend(loc="upper right")

    plt.savefig(os.path.join(scan_dir,"spectra.png"),dpi=1000)
    plt.show()
