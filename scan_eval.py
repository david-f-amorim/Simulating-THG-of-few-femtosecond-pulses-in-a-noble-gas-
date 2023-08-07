import numpy as np
import matplotlib.pyplot as plt
import os  

"""
 scan_eval.py: 

        Processing and visualisation of pressure scan output. 
"""

# ---------- EXECUTION SETTINGS --------------------------------------

sup_dir    = "scans"
out_dir    = "scan_analysis"


single_dir = "scans\scan_150.0mW_Ne_0.0rad_f_ion" 
single = False 
n = 15

# ---------- FILE EXTRACTION --------------------------------------

# get paths to all files in super directory with different 
# beam powers and all other parameters equal
def power_comp(sup_dir, gas, phi,ion, kerr):

    path_arr =np.array([])
    beam_en_arr = np.array([])

    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)):
            
            gas0, phi0, beam_en, ion0, kerr0 = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (phi0 == phi) & (kerr0 == kerr) & (ion0 == ion) 

            if cond:
                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                beam_en_arr = np.append(beam_en_arr, beam_en)

    
    path_arr = path_arr[beam_en_arr.argsort()] 
    beam_en_arr = beam_en_arr[beam_en_arr.argsort()]
             

    return path_arr, beam_en_arr 

# get paths to all files in super directory with different 
# CEO phases and all other parameters equal
def phi_comp(sup_dir, gas, beam_en, ion,kerr):

    path_arr =np.array([])
    phi_arr = np.array([])

    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)):
            
            gas0, phi, beam_en0, ion0, kerr0 = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (beam_en0 == beam_en) & (kerr0 == kerr) & (ion0 == ion) 

            if cond:
                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                phi_arr = np.append(phi_arr, phi)

    
    path_arr = path_arr[phi_arr.argsort()] 
    phi_arr = phi_arr[phi_arr.argsort()]

    return path_arr, phi_arr


# get paths to all files in super directory with different 
# beam energies and all other parameters equal
def gas_comp(sup_dir,phi,beam_en,ion, kerr):

    path_arr = []

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            _,phi0, beam_en0, ion0, kerr0 = get_params(os.path.join(sup_dir,dir))

            cond = (beam_en0==beam_en) & (phi0 == phi) & (kerr0 == kerr) & (ion0 == ion) 

            if cond:
                path_arr.append(os.path.join(sup_dir,dir))

    return path_arr 

# get paths to all files in super directory with different 
# ionisation and all other parameters [apart from power!] equal
# (only get path to ion file if corresponding non-ion file exists)
def ion_comp(sup_dir, gas, phi, kerr):

    path_arr_no_ion = np.array([])
    beam_en_arr = np.array([])
    path_arr_ion = np.array([])

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            gas0, phi0, beam_en, ion, kerr0 = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (phi0 == phi) & (kerr0 == kerr) & (ion=="false")

            if cond:
                path_arr_no_ion = np.append(path_arr_no_ion,os.path.join(sup_dir,dir))
                beam_en_arr = np.append(beam_en_arr,beam_en)

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            gas0, phi0, beam_en, ion, kerr0 = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (phi0 == phi) & (kerr0 == kerr) & (ion=="true") & (beam_en in beam_en_arr)

            if cond:
                path_arr_ion =np.append(path_arr_ion,os.path.join(sup_dir,dir))
           

    path_arr_no_ion = path_arr_no_ion[beam_en_arr.argsort()] 
    path_arr_ion = path_arr_ion[beam_en_arr.argsort()]
    beam_en_arr = beam_en_arr[beam_en_arr.argsort()]

    return path_arr_no_ion, path_arr_ion, beam_en_arr 

# ---------- DATA EXTRACTION --------------------------------------

# extract parameter values of a scan directory at 'path' 
def get_params(path):

    params = np.loadtxt(os.path.join(path,"params.txt"),dtype="str", delimiter="=", comments="#")

    gas  = params[0,1]
    phi = float(params[6,1])
    beam_en = float(params[7,1])
    kerr = params[9,1]
    ion = params[10,1]

    return gas[1:], phi, beam_en, ion[1:], kerr[1:]  

# extract pressures, energies, and efficiencies of a scan 
# directory at 'path'
def get_PrEnEf(path):

    arr    = np.loadtxt(os.path.join(path,"energy_efficiency.txt"))
    arr    = arr[arr[:, 0].argsort()]
    p_arr  = arr[:,0] 
    en_arr = arr[:,1]
    ef_arr = arr[:,2]

    return p_arr, en_arr, ef_arr

# extract peak pressure, peak energy, and peak efficiency of 
# a scan directory at 'path'
def get_peaks(path):

    p_arr, en_arr, ef_arr = get_PrEnEf(path)

    peak_en = np.max(en_arr)
    peak_ef = np.max(ef_arr)
    peak_p =  p_arr[np.where(en_arr == peak_en )][0]

    return peak_p, peak_en, peak_ef

# extract spectra from scan directory at 'path'
def get_spectra(path):

    p_arr, _, _ = get_PrEnEf(path)
    N = len(p_arr)
    data = np.empty(shape=(N,3), dtype="object")

    for i in np.arange(N):

        tmp_arr = np.loadtxt(os.path.join(path,str(p_arr[i])+" bar.dat"))
        lam = tmp_arr[:,0]
        I = tmp_arr[:,1]

        data[i,0] = p_arr[i]
        data[i,1] = lam 
        data[i,2] = I 

    return data 

# reduce number of spectra to n:
def reduce_spectra(path,n):
    data = get_spectra(path)
    p_arr = data[:,0]

    N = len(data[:,0])

    if N <= n:
        return data 

    data_cut = np.empty(shape=(n,3), dtype="object")
    p_cut = np.empty(n)
    k = np.max(p_arr)/n

    for i in np.arange(n):

        data_cut[i,:] = data[np.argmin(np.abs(p_arr-i*k)),:]
        p_cut[i] = p_arr[np.argmin(np.abs(p_arr-i*k))]

    return data_cut, p_cut     


# ---------- VISUALISATION --------------------------------------

# plot UV energy, THG conversion efficiency and UV spectra for a 
# single scan 
def plot_single(single_dir, n=15):

    p_arr, en_arr, ef_arr = get_PrEnEf(single_dir)
    p_peak, en_peak, ef_peak = get_peaks(single_dir)
    gas, phi, beam_en, ion, kerr = get_params(single_dir)

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; beam energy: {0}mW; CEO phase: {1:.2f}rad; ".format(beam_en*1e6, phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, en_arr*1e9, color="blue")
    plt.scatter(p_arr, en_arr*1e9, color="blue", label="Peak: {0:.1f}nJ at {1}bar".format(en_peak*1e9, p_peak))
    plt.legend()

    plt.savefig(os.path.join(single_dir,"energies.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; beam energy: {0}mW; CEO phase: {1:.2f}rad; ".format(beam_en*1e6, phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, ef_arr*1e2, color="blue")
    plt.scatter(p_arr, ef_arr*1e2, color="blue", label="Peak: {0:.2f}% at {1}bar".format(ef_peak*1e2, p_peak))
    plt.legend()

    plt.savefig(os.path.join(single_dir,"efficiencies.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Simulated UV spectra", fontsize=16)
    plt.title("Gas: "+gas+"; beam energy: {0}mW; CEO phase: {1:.2f}rad; ".format(beam_en*1e6, phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10, pad=18)
    plt.ylabel("Intensity (arb. units)")
    plt.xlabel("Wavelength (nm)")

    N = len(p_arr)
    cmap = plt.get_cmap("viridis")
    data, p_cut = reduce_spectra(single_dir, n)
    cidx = p_cut / np.max(p_cut) 
    
    for i in np.arange(n):
        plt.plot(data[i,1]*1e9, data[i,2], color=cmap(cidx[i]), label="{0}bar".format(data[i,0]))
    plt.legend(loc="upper right")

    plt.savefig(os.path.join(single_dir,"spectra.png"),dpi=1000)
    plt.show()

# plot UV energy and THG conversion efficiency for different beam power
def plot_beamP_scan(sup_dir, gas, phi, ion, kerr):

    # create output directory
    if ion=="true":
        ion_string ="ion"
    else:    
        ion_string="no-ion"    

    out_path = os.path.join(out_dir, "beamP_scan_"+gas+"_"+str(phi)+"rad_"+kerr+"_"+ion_string)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr, beam_en_arr = power_comp(sup_dir, gas, phi,ion, kerr)

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = beam_en_arr / np.max(beam_en_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _ = get_PrEnEf(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW".format(beam_en_arr[i]*1e6))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr = get_PrEnEf(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW".format(beam_en_arr[i]*1e6))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    plt.show()

    # PLOTS 3-5: peak investigation 
    peak_arr = np.empty((len(path_arr),4))

    for i in np.arange(len(path_arr)):
        peak_arr[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak = get_peaks(path_arr[i])
        peak_arr[i, 1] = p_peak 
        peak_arr[i,2] = en_peak 
        peak_arr[i,3] = ef_peak 

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Peak UV energy (nJ)")
    plt.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")
    plt.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Efficiency (%)")
    plt.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
    plt.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Saturation pressure (bar)")
    plt.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
    plt.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")

    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    plt.show()

# plot UV energy and THG conversion efficiency for different beam powers
# with and without ionisation
def plot_beamP_ion_scan(sup_dir, gas, phi, kerr):

    # create output directory    

    out_path = os.path.join(out_dir, "beamP_ion_scan_"+gas+"_"+str(phi)+"rad_"+kerr)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr_no_ion, path_arr_ion, beam_en_arr = ion_comp(sup_dir, gas, phi, kerr)

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = beam_en_arr / np.max(beam_en_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, en_arr, _ = get_PrEnEf(path_arr_ion[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW (ion.)".format(beam_en_arr[i]*1e6))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

        p_arr, en_arr, _ = get_PrEnEf(path_arr_no_ion[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW (no ion.)".format(beam_en_arr[i]*1e6), marker="+")
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]), ls="--")

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, _, ef_arr = get_PrEnEf(path_arr_ion[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW (ion.)".format(beam_en_arr[i]*1e6))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

        p_arr, _, ef_arr = get_PrEnEf(path_arr_no_ion[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW (no ion.)".format(beam_en_arr[i]*1e6), marker="+")
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    plt.show()

    # PLOTS 3-5: peak investigation 
    peak_arr_ion = np.empty((len(path_arr_ion),4))
    peak_arr_no_ion = np.empty((len(path_arr_ion),4))

    for i in np.arange(len(path_arr_ion)):
        peak_arr_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak = get_peaks(path_arr_ion[i])
        peak_arr_ion[i, 1] = p_peak 
        peak_arr_ion[i,2] = en_peak 
        peak_arr_ion[i,3] = ef_peak 

        peak_arr_no_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak = get_peaks(path_arr_no_ion[i])
        peak_arr_no_ion[i, 1] = p_peak 
        peak_arr_no_ion[i,2] = en_peak 
        peak_arr_no_ion[i,3] = ef_peak 

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Peak UV energy (nJ)")
    plt.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue", label="ion.")
    plt.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue")
    plt.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red", label="no ion.")
    plt.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red")

    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Efficiency (%)")
    plt.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue", label="ion.")
    plt.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue")
    plt.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red", label="no ion.")
    plt.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red")
    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; CEO phase: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    plt.xlabel("Beam power (mW)")
    plt.ylabel("Saturation pressure (bar)")
    plt.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue", label="ion.")
    plt.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue")
    plt.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red", label="no ion.")
    plt.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red")

    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    plt.show()    

# plot UV energy and THG conversion efficiency for different beam power
def plot_phi_scan(sup_dir, gas, beam_en, ion, kerr):

    # create output directory
    if ion=="true":
        ion_string ="ion"
    else:    
        ion_string="no-ion"    

    out_path = os.path.join(out_dir, "phi_scan_"+gas+"_{0}mW_".format(beam_en*1e6)+kerr+"_"+ion_string)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr, phi_arr = phi_comp(sup_dir, gas, beam_en,ion, kerr)

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = phi_arr / np.max(phi_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW; ".format(beam_en*1e6)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _ = get_PrEnEf(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW; ".format(beam_en*1e6)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr = get_PrEnEf(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    plt.show()

    # PLOTS 3-5: peak investigation 
    peak_arr = np.empty((len(path_arr),4))

    for i in np.arange(len(path_arr)):
        peak_arr[i,0] = phi_arr[i]
        p_peak, en_peak, ef_peak = get_peaks(path_arr[i])
        peak_arr[i, 1] = p_peak 
        peak_arr[i,2] = en_peak 
        peak_arr[i,3] = ef_peak 

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW; ".format(beam_en*1e6)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    plt.xlabel("CEO phase (rad)")
    plt.ylabel("Peak UV energy (nJ)")
    plt.scatter(peak_arr[:,0], peak_arr[:,2]*1e9, color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,2]*1e9, color="blue")

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW; ".format(beam_en*1e6)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    plt.xlabel("CEO phase (rad)")
    plt.ylabel("Efficiency (%)")
    plt.scatter(peak_arr[:,0], peak_arr[:,3]*1e2, color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,3]*1e2, color="blue")

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW; ".format(beam_en*1e6)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    plt.xlabel("CEO phase (rad)")
    plt.ylabel("Saturation pressure (bar)")
    plt.scatter(peak_arr[:,0], peak_arr[:,1], color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,1], color="blue")

    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    plt.show()

# ---------- EXEC --------------------------------------
if single:
    plot_single(single_dir, n)
else:
    plot_beamP_scan(sup_dir, "Ar", 0.0, "true", "f") 
    plot_beamP_scan(sup_dir, "Ar", 0.0, "false", "f") 
    plot_beamP_ion_scan(sup_dir, "Ar", 0.0, "f") 
    plot_phi_scan(sup_dir,"Ar", 150e-6, "true", "f")

    plot_beamP_scan(sup_dir, "Ne", 0.0, "true", "f") 
    plot_beamP_scan(sup_dir, "Ne", 0.0, "false", "f") 
    plot_beamP_ion_scan(sup_dir, "Ne", 0.0, "f") 
    


