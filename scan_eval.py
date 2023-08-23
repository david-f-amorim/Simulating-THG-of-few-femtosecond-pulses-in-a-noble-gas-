import numpy as np
import matplotlib.pyplot as plt
import os  
import json
from datetime import DateTime
# ---------- QUICK SETTINGS --------------------------------------

single = False  # if True: process output from a single pressure scan; if False: process output from a set of different pressure scans ("multi-scan")

n    = 15   # maximum number of overlayed spectra in one plot (to avoid clutter)
show = True # if True: open plots in matplotlib GUI; if False: don't show plots (directly write to file)  

# ---------- MULTI-SCAN SETTINGS --------------------------------------



# ---------- INPUT/OUTPUT HANDLING --------------------------------------

single_dir = "scan_new_prelims\\scan_300.0mW_Ar_0.0rad_f_ion_grad" # path to pressure scan directory to process if single=True (Note: output will be written to same directory)

sup_dir = "parameter_scans\\gas_scans" # if single==False, path to super directory containing the various pressure scan directories 
out_dir = "scan_analysis\\gas_scans"   # if single==False, path to output directory for the produced plots 

# ---------- MULTI-SCANS: FILE EXTRACTION --------------------------------------

# get paths to all pressure scan files in sup_dir with different 
# beam powers and all other parameters equal
#    NOTE: if only sup_dir is specified, this assumes that all 
#          pressure scans in sup_dir have all parameters equal 
#          except beam power; if this is not the case, provide 
#          a function argument to filter for relevant parameters
#          using **kwargs (accepted kwargs: all parameters 
#          returned by get_params)

def power_comparison_get(sup_dir, **kwargs ):

    path_arr =np.array([])
    IR_energy_arr = np.array([])
    IR_int_arr = np.array([])

    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)):
            
            params_arr, params_dict = get_params(os.path.join(sup_dir,dir))

            cond = 1

            for key, value in kwargs.items():
                if (key in params_dict) and (value==params_dict.get(key)):
                    cond *=1
                else:
                    cond *=0   

            if cond: 
                IR_energy = params_arr[6]
                IR_int    = params_arr[12]

                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                IR_energy_arr = np.append(IR_energy_arr, IR_energy)
                IR_int_arr = np.append(IR_int_arr, IR_int)

    
    path_arr = path_arr[IR_energy_arr.argsort()] 
    IR_int_arr = IR_int_arr[IR_energy_arr.argsort()]
    IR_energy_arr = IR_energy_arr[IR_energy_arr.argsort()]
             

    return path_arr, IR_energy_arr, IR_int_arr

# get paths to all files in super directory with different 
# CEPs and all other parameters equal
def phi_comp(sup_dir, gas, beam_en, ion, dens_mod):

    path_arr =np.array([])
    phi_arr = np.array([])

    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
            
            gas0, phi, beam_en0, ion0, kerr0, _, _, I = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (beam_en0 == beam_en) & (kerr0 == kerr) & (ion0 == ion) 

            if cond:
                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                phi_arr = np.append(phi_arr, phi)
    
    path_arr = path_arr[phi_arr.argsort()] 
    phi_arr = phi_arr[phi_arr.argsort()]

    return path_arr, phi_arr, I

# get paths to all files in super directory with different 
# beam energies and all other parameters equal
def gas_comp_singleP(sup_dir,beam_en,phi, ion, dens_mod):

    path_arr = []
    gas_arr = []

    for dir in os.listdir(sup_dir) :
        if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
            gas,_, beam_en0, _, _, _, _, I = get_params(os.path.join(sup_dir,dir))

            cond = (beam_en0==beam_en*1e-3)

            if cond:
                path_arr.append(os.path.join(sup_dir,dir))
                gas_arr.append(gas)

    return path_arr, gas_arr 

# get paths to all files in super directory with different 
# ionisation and all other parameters [apart from power!] equal
# (only get path to ion file if corresponding non-ion file exists)
def ion_comp(sup_dir, gas, phi, dens_mod):

    path_arr_no_ion = np.array([])
    beam_en_arr = np.array([])
    I_arr = np.array([])
    path_arr_ion = np.array([])

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
            gas0, phi0, beam_en, ion, kerr0, w, tau, I = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (phi0 == phi) & (kerr0 == kerr) & (ion=="false")

            if cond:
                path_arr_no_ion = np.append(path_arr_no_ion,os.path.join(sup_dir,dir))
                beam_en_arr = np.append(beam_en_arr,beam_en)
                I_arr = np.append(I_arr, I)

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
            gas0, phi0, beam_en, ion, kerr0, w, tau, I = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas) & (phi0 == phi) & (kerr0 == kerr) & (ion=="true") & (beam_en in beam_en_arr)

            if cond:
                path_arr_ion =np.append(path_arr_ion,os.path.join(sup_dir,dir))
           

    path_arr_no_ion = path_arr_no_ion[beam_en_arr.argsort()] 
    path_arr_ion = path_arr_ion[beam_en_arr.argsort()]
    I_arr = I_arr[beam_en_arr.argsort()]
    beam_en_arr = beam_en_arr[beam_en_arr.argsort()]

    return path_arr_no_ion, path_arr_ion, beam_en_arr, I_arr, w, tau 

# for a single gas and a single power [in mW] show both dens models
def singlePg_model_comp(sup_dir, gas, ion, phi, beam_p):

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            gas0, _, beam_en, _, _, w, tau, I = get_params(os.path.join(sup_dir,dir))

            cond_coms = (gas0==gas) & (beam_en*1e3==beam_p) & (dir[-4:]=="coms")
            cond_grad = (gas0==gas) & (beam_en*1e3==beam_p) & (dir[-4:]=="grad")

            if cond_coms:
                path_coms =os.path.join(sup_dir,dir)
            elif cond_grad:
                path_grad =os.path.join(sup_dir,dir)
             

    return path_coms, path_grad, beam_en, w, tau      

# compare power scans of single gas with different dens models 
def multiP_singleg_model_comp(sup_dir, gas):

    path_arr_coms = np.array([])
    beam_en_arr = np.array([])
    I_arr = np.array([])
    path_arr_grad = np.array([])

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            gas0, _, beam_en, _, _, w, tau, I = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas)  & (dir[-4:]=="coms") 

            if cond:
                path_arr_coms = np.append(path_arr_coms,os.path.join(sup_dir,dir))
                beam_en_arr = np.append(beam_en_arr,beam_en)
                I_arr = np.append(I_arr,I)
                

    for dir in os.listdir(sup_dir):
        if os.path.isdir(os.path.join(sup_dir,dir)):
            gas0, _, beam_en, _, _, w, tau, I = get_params(os.path.join(sup_dir,dir))

            cond = (gas0==gas)  & (dir[-4:]=="grad") & (beam_en in beam_en_arr)

            if cond:
                path_arr_grad =np.append(path_arr_grad,os.path.join(sup_dir,dir))
                
           

    path_arr_coms = path_arr_coms[beam_en_arr.argsort()] 
    path_arr_grad = path_arr_grad[beam_en_arr.argsort()]
    I_arr = I_arr[beam_en_arr.argsort()]
    beam_en_arr = beam_en_arr[beam_en_arr.argsort()]

    return path_arr_grad, path_arr_coms, beam_en_arr, I_arr 
    
# get paths to all files in super directory with different 
# beam powers and all other parameters equal
def gas_comp_multiP(sup_dir,dens_mod, excluded_gases):

    gas_arr = np.array([])
  
    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
            
            gas, _, beam_en,_, _, w0, tau, I = get_params(os.path.join(sup_dir,dir))

            if not (gas in gas_arr) and not (gas in excluded_gases):
                gas_arr = np.append(gas_arr, gas)

    data = np.empty(shape=(len(gas_arr), 3), dtype="object")
    
    for i in np.arange(len(gas_arr)):
        path_arr =np.array([])
        beam_en_arr = np.array([])
        for dir in os.listdir(sup_dir):

            if os.path.isdir(os.path.join(sup_dir,dir)) and (dir[-4:]==dens_mod):
                
                gas0, _, beam_en,_, _, w0, tau, I = get_params(os.path.join(sup_dir,dir))

                cond = (gas0==gas_arr[i]) 

                if cond:
                    path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                    beam_en_arr = np.append(beam_en_arr, beam_en)

        path_arr = path_arr[beam_en_arr.argsort()] 
        beam_en_arr = beam_en_arr[beam_en_arr.argsort()]

        data[i, 0]=gas_arr[i]
        data[i, 1]= path_arr 
        data[i, 2]= beam_en_arr
    
    return data, w0, tau     

# ---------- FUNCTIONS THAT EXTRACT DATA FROM A SINGLE PRESSURE SCAN DIRECTORY --------------------------------------

# extract parameter values of the pressure scan directory at 'path' and returns array and dict of parameters
def get_params(path):

    params = np.loadtxt(os.path.join(path,"params.txt"),dtype="str", delimiter="=", comments="#")

    gas  = params[0,1][1:]
    dens_mod =params[2,1]
    tau = float(params[3,1])
    lam0=float(params[4,1])
    w0 = float(params[5,1])
    CEP = float(params[6,1])
    IR_energy = float(params[7,1])
    ion = params[10,1][1:]
    propz=float(params[11,1])
    GVD = float(params[12,1])
    thickness=float(params[13,1])
    ion_mod=params[14,1][1:]
    
    # also calculate 1/e^2 IR intensity in W/cm^2:
    IR_int = IR_energy  / (np.pi * (w0*1e2)**2 * tau)

    params_arr = np.array([gas, dens_mod, tau, lam0, w0, CEP, IR_energy, ion, propz, GVD, thickness, ion_mod, IR_int]) 
    params_dict= {"gas": gas, "dens_mod": dens_mod, "tau": tau, "lam0": lam0, "w0": w0, "CEP": CEP, "IR_energy": IR_energy, "ion": ion, "propz": propz, "GVD": GVD, "thickness": thickness, "ion_mod": ion_mod, "IR_int": IR_int}

    return params_arr, params_dict

# extract pressures, energies, efficiencies, pulse duraton, and peak position of the scan  directory at 'path'; 
# returns different arrays
def get_data(path):

    arr    = np.loadtxt(os.path.join(path,"energy_efficiency_time_zpeak.txt"))
    arr    = arr[arr[:, 0].argsort()]
    p_arr  = arr[:,0] 
    en_arr = arr[:,1]
    ef_arr = arr[:,2]
    tau_arr= arr[:,3]
    zpeak_arr= arr[:,4]

    return p_arr, en_arr, ef_arr, tau_arr, zpeak_arr 

# extract peak energy, peak efficiency, minimum pulse duration, pressure at peak energy, and pressure at minimum pulse duration of the pressure scan directory at 'path';
# returns array 
def get_peaks(path):

    p_arr, en_arr, ef_arr, tau_arr, _ = get_data(path)

    peak_en = np.max(en_arr)
    peak_ef = np.max(ef_arr)
    min_tau = np.min(tau_arr)
    min_tau_pres =p_arr[np.where(tau_arr == min_tau )][0]
    peak_en_pres = p_arr[np.where(en_arr == peak_en )][0]

    return np.array([peak_en, peak_ef, min_tau, peak_en_pres, min_tau_pres])

# extract spectra from pressure scan directory at 'path' and reduce number of spectra to n ; returns arrays of spectra and corresponding pressure values
def get_spectra(path, n):

    p_arr, _, _, _, _ = get_data(path)
    N = len(p_arr)
    data = np.empty(shape=(N,3), dtype="object")

    for i in np.arange(N):

        tmp_arr = np.loadtxt(os.path.join(path,str(p_arr[i])+" bar.dat"))
        lam = tmp_arr[:,0]
        I = tmp_arr[:,1]

        data[i,0] = p_arr[i]
        data[i,1] = lam 
        data[i,2] = I 

    if N <= n:
        return data, p_arr 
    
    else:
        data_cut = np.empty(shape=(n,3), dtype="object")
        p_cut = np.empty(n)
        k = np.max(p_arr)/n

        for i in np.arange(n):

            data_cut[i,:] = data[np.argmin(np.abs(p_arr-i*k)),:]
            p_cut[i] = p_arr[np.argmin(np.abs(p_arr-i*k))]

        return data_cut, p_cut

# ---------- VISUALISATION --------------------------------------

# define set of colours 
colour_cycle = ["grey", "black", "red", "blue", "purple", "green", "cyan", "orange", "deeppink"]

# plot UV energy, THG conversion efficiency and UV spectra for a  single scan 
def plot_single(single_dir, n=15):

    p_arr, en_arr, ef_arr, tau_arr, zpeak_arr = get_data(single_dir)
    p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(single_dir)
    gas, phi, beam_en, ion, kerr, _, _, I = get_params(single_dir)

    dens_mod = single_dir[-4:]

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10, pad=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, en_arr*1e9, color="blue")
    plt.scatter(p_arr, en_arr*1e9, color="blue", label="Peak: {0:.1f}nJ at {1}bar".format(en_peak*1e9, p_peak))
    plt.legend(loc="upper left")

    plt.savefig(os.path.join(single_dir,"energies.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10, pad=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, ef_arr*1e2, color="blue")
    plt.scatter(p_arr, ef_arr*1e2, color="blue", label="Peak: {0:.2f}% at {1}bar".format(ef_peak*1e2, p_peak))
    plt.legend()

    plt.savefig(os.path.join(single_dir,"efficiencies.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated pulse durations", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10, pad=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, tau_arr*1e15, color="blue")
    plt.scatter(p_arr, tau_arr*1e15, color="blue", label="Minimum: {0:.2f}fs at {1}bar".format(tau_peak*1e15, min_tau_p))
    plt.legend()

    plt.savefig(os.path.join(single_dir,"durations.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10, pad=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, zpeak_arr*1e3, color="blue")
    plt.scatter(p_arr, zpeak_arr*1e3, color="blue", label="Maximum: {0:.2f}mm at {1}bar".format(np.max(zpeak_arr)*1e3,p_arr[np.where(zpeak_arr == np.max(zpeak_arr) )][0]))
    plt.legend()

    plt.savefig(os.path.join(single_dir,"z_peak.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.82)
    plt.suptitle("Simulated UV spectra", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I(single_dir)*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10, pad=18)
    plt.ylabel("Intensity (arb. units)")
    plt.xlabel("Wavelength (nm)")

    N = len(p_arr)
    cmap = plt.get_cmap("viridis")
    data, p_cut = get_spectra(single_dir, n)
    cidx = p_cut / np.max(p_cut) 
    
    if N<=n:
        n=N

    for i in np.arange(n):
        plt.plot(data[i,1]*1e9, data[i,2], color=cmap(cidx[i]), label="{0}bar".format(data[i,0]))
    plt.legend(loc="upper right")

    plt.savefig(os.path.join(single_dir,"spectra.png"),dpi=1000)
    if show: plt.show()

# compare single gas, single beam power, two dens models
def plot_singlePg_model_comp(sup_dir, gas, beam_p):

    path_coms, path_grad, beam_en, w, tau = singlePg_model_comp(sup_dir, gas, beam_p)

    p_arr_coms, en_arr_coms, ef_arr_coms, tau_arr_coms, zpeak_arr_coms = get_data(path_coms)
    p_peak_coms, en_peak_coms, ef_peak_coms, tau_peak_coms, min_tau_p_coms = get_peaks(path_coms)
    p_arr_grad, en_arr_grad, ef_arr_grad, tau_arr_grad, zpeak_arr_grad = get_data(path_grad)
    p_peak_grad, en_peak_grad, ef_peak_grad, tau_peak_grad, min_tau_p_grad = get_peaks(path_grad)

    gas, phi, beam_en, ion, kerr, _, _, I = get_params(path_coms)


    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I(single_dir)*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion, fontsize=10, pad=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr_coms, en_arr_coms*1e9, color="blue")
    plt.scatter(p_arr_coms, en_arr_coms*1e9, color="blue", label="COMSOL (peak: {0:.1f}nJ at {1}bar)".format(en_peak_coms*1e9, p_peak_coms))
    plt.plot(p_arr_grad, en_arr_grad*1e9, color="red")
    plt.scatter(p_arr_grad, en_arr_grad*1e9, color="red", label="gradient (peak: {0:.1f}nJ at {1}bar)".format(en_peak_grad*1e9, p_peak_grad))
    plt.legend()

    plt.savefig(os.path.join(path_coms,"energies.png"),dpi=1000)
    plt.savefig(os.path.join(path_grad,"energies.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion, fontsize=10, pad=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr_coms, ef_arr_coms*1e2, color="blue")
    plt.scatter(p_arr_coms, ef_arr_coms*1e2, color="blue", label="COMSOL (peak: {0:.2f}% at {1}bar)".format(ef_peak_coms*1e2, p_peak_coms))
    plt.plot(p_arr_grad, ef_arr_grad*1e2, color="red")
    plt.scatter(p_arr_grad, ef_arr_grad*1e2, color="red", label="gradient (peak: {0:.2f}% at {1}bar)".format(ef_peak_grad*1e2, p_peak_grad))
    plt.legend()

    plt.savefig(os.path.join(path_coms,"efficiencies.png"),dpi=1000)
    plt.savefig(os.path.join(path_grad,"efficiencies.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.84)
    plt.suptitle("Simulated pulse durations", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {2:.1f}PW/cm^2); CEP: {1:.2f}rad; ".format(beam_en*1e6, phi, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion, fontsize=10, pad=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr_coms, tau_arr_coms*1e15, color="blue")
    plt.scatter(p_arr_coms, tau_arr_coms*1e15, color="blue", label="COMSOL (minimum: {0:.2f}fs at {1}bar)".format(tau_peak_coms*1e15, min_tau_p_coms))
    plt.plot(p_arr_grad, tau_arr_grad*1e15, color="red")
    plt.scatter(p_arr_grad, tau_arr_grad*1e15, color="red", label="gradient (minimum: {0:.2f}fs at {1}bar)".format(tau_peak_grad*1e15, min_tau_p_grad))
    plt.legend()

    plt.savefig(os.path.join(path_coms,"durations.png"),dpi=1000)
    plt.savefig(os.path.join(path_grad,"durations.png"),dpi=1000)
    if show: plt.show()

# peak investigations for multi power, single gas, dens model comparison 
def plot_multiP_singleg_model_comp(sup_dir, gas):

    out_path = os.path.join(out_dir, "beamP_scan_"+gas+"_dens_comp")
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr_grad, path_arr_coms, beam_en_arr, I_arr =multiP_singleg_model_comp(sup_dir, gas)

    # get trivia
    _, phi, _, ion, kerr, w0, tau, I = get_params(path_arr_coms[0])

     # auxiliary functions (beam_en in mW, I in PW/cm^2):
    def p2i(beam_en):
        return beam_en *1e-6 / (np.pi * (w0*1e2)**2 * tau ) * 1e-15

    def i2p(I):
        return I * (np.pi * (w0*1e2)**2 * tau ) *1e15 *1e6    
    
    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = beam_en_arr / np.max(beam_en_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_coms)):
        p_arr, en_arr, _, _, _ = get_data(path_arr_coms[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (COMSOL)".format(beam_en_arr[i]*1e6, 1e-15*I_arr[i]))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

        p_arr, en_arr, _, _, _ = get_data(path_arr_grad[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (gradient)".format(beam_en_arr[i]*1e6, 1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]), ls="--")

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_coms)):
        p_arr, _, ef_arr,_, _ = get_data(path_arr_coms[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2]  (COMSOL)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

        p_arr, _, ef_arr, _, _ = get_data(path_arr_grad[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (gradient)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV pulse duration", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_coms)):
        p_arr, _, _,tau_arr, _ = get_data(path_arr_coms[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2]  (COMSOL)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]))

        p_arr, _, _, tau_arr, _ = get_data(path_arr_grad[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (gradient)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

   
    # PLOTS 5-9: peak investigation 
    peak_arr_ion = np.empty((len(path_arr_coms),6))
    peak_arr_no_ion = np.empty((len(path_arr_coms),6))

    for i in np.arange(len(path_arr_coms)):
        peak_arr_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr_coms[i])
        peak_arr_ion[i, 1] = p_peak 
        peak_arr_ion[i,2] = en_peak 
        peak_arr_ion[i,3] = ef_peak 
        peak_arr_ion[i,4] = tau_peak
        peak_arr_ion[i,5] = min_tau_p

        peak_arr_no_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr_grad[i])
        peak_arr_no_ion[i, 1] = p_peak 
        peak_arr_no_ion[i,2] = en_peak 
        peak_arr_no_ion[i,3] = ef_peak 
        peak_arr_no_ion[i,4] = tau_peak
        peak_arr_no_ion[i,5] = min_tau_p

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue", label="COMSOL")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red", label="gradient")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue", label="COMSOL")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red", label="gradient")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue", label="COMSOL")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red", label="gradient")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()  

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated UV pulse durations ", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Mininum pulse duration (fs)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,4]*1e15, color="blue", label="COMSOL")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,4]*1e15, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,4]*1e15, color="red", label="gradient")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,4]*1e15, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    plt.savefig(os.path.join(out_path,"min_tau_vs_beam_power.png"),dpi=1000)
    if show: plt.show()   

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated minimum pulse duration pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Mininum pulse duration pressure (bar)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,5], color="blue", label="COMSOL")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,5], color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,5], color="red", label="gradient")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,5], color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

# plot UV energy and THG conversion efficiency for different beam power
def power_comparison_plot(sup_dir, **kwargs):

    # get data and params
    path_arr, IR_energy_arr, IR_int_arr = power_comparison_get(sup_dir, **kwargs)

    params, params_dict = get_params(path_arr[0])
    w0, tau = params[4], params[2]

    # set title string 
    title_str = ""
    for key, val in kwargs.items():
        title_str += key+": "+val+"; "
  
    # set output directory 
    out_path = os.path.join(out_dir, "beam_power_comparison_{}") and "add date time here!!"
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # auxiliary functions (beam_en in mW, I in PW/cm^2):
    def p2i(beam_en):
        return beam_en *1e-6 / (np.pi * (w0*1e2)**2 * tau ) * 1e-15

    def i2p(I):
        return I * (np.pi * (w0*1e2)**2 * tau ) *1e15 *1e6
    
    # set title string

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = beam_en_arr / np.max(beam_en_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title(title_str, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _,_,_ = get_data(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW (I={1:.1f}PW/cm^2)".format(beam_en_arr[i]*1e6, 1e-15*I_arr[i]))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title(title_str, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr, _, _ = get_data(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW (I={1:.1f}PW/cm^2)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: pulse duration vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV pulse durations", fontsize=16)
    plt.title(title_str, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, tau_arr,_ = get_data(path_arr[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0}mW (I={1:.1f}PW/cm^2)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

    # PLOT 4: z_peak vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title(title_str, fontsize=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, _,z_peak_arr = get_data(path_arr[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]), label="{0}mW (I={1:.1f}PW/cm^2)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_positions.png"),dpi=1000)
    if show: plt.show()

    # PLOTS 4-7: peak investigation 
    peak_arr = np.empty((len(path_arr),6))

    for i in np.arange(len(path_arr)):
        peak_arr[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr[i])
        peak_arr[i, 1] = p_peak 
        peak_arr[i,2] = en_peak 
        peak_arr[i,3] = ef_peak 
        peak_arr[i,4] = tau_peak 
        peak_arr[i,5] = min_tau_p

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")
    ax.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")
    ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
    ax.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
    ax.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    plt.suptitle("Simulated UV pulse durations", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration (fs)")
    ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue")
    ax.plot(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"tau_min_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    plt.suptitle("Simulated minimal pulse duration pressure", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration pressure (bar)")
    ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue")
    ax.plot(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    # save kwargs and sup_dir to file
    params_dict["sup_dir"] = sup_dir

    with open(os.path.join(out_path,"kwargs.txt"), "w") as file:
        json.dump(params_dict, file)



# plot UV energy and THG conversion efficiency for different beam powers
# with and without ionisation
def plot_beamP_ion_scan(sup_dir, gas, phi, kerr, dens_mod):

    # create output directory    

    out_path = os.path.join(out_dir, "beamP_ion_scan_"+gas+"_"+str(phi)+"rad_"+kerr+"_"+dens_mod)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr_no_ion, path_arr_ion, beam_en_arr, I_arr, w0, tau = ion_comp(sup_dir, gas, phi, kerr, dens_mod)

     # auxiliary functions (beam_en in mW, I in PW/cm^2):
    def p2i(beam_en):
        return beam_en *1e-6 / (np.pi * (w0*1e2)**2 * tau ) * 1e-15

    def i2p(I):
        return I * (np.pi * (w0*1e2)**2 * tau ) *1e15 *1e6

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = beam_en_arr / np.max(beam_en_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, en_arr, _, _, _ = get_data(path_arr_ion[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (ion.)".format(beam_en_arr[i]*1e6, 1e-15*I_arr[i]))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

        p_arr, en_arr, _, _, _ = get_data(path_arr_no_ion[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (no ion.)".format(beam_en_arr[i]*1e6, 1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]), ls="--")

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, _, ef_arr,_, _ = get_data(path_arr_ion[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2]  (ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

        p_arr, _, ef_arr, _, _ = get_data(path_arr_no_ion[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (no ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV pulse duration", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, _, _,tau_arr, _ = get_data(path_arr_ion[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2]  (ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]))

        p_arr, _, _, tau_arr, _ = get_data(path_arr_no_ion[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (no ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr_ion)):
        p_arr, _, _,_, z_peak_arr = get_data(path_arr_ion[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2]  (ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]))
        plt.plot(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]))

        p_arr, _, _, _, z_peak_arr = get_data(path_arr_no_ion[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]), label="{0}mW [I={1:.1f} PW/cm^2] (no ion.)".format(beam_en_arr[i]*1e6,1e-15*I_arr[i]), marker="+")
        plt.plot(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]), ls='--')

    plt.legend()
    plt.savefig(os.path.join(out_path,"z_peak.png"),dpi=1000)
    if show: plt.show()

    # PLOTS 5-9: peak investigation 
    peak_arr_ion = np.empty((len(path_arr_ion),6))
    peak_arr_no_ion = np.empty((len(path_arr_ion),6))

    for i in np.arange(len(path_arr_ion)):
        peak_arr_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr_ion[i])
        peak_arr_ion[i, 1] = p_peak 
        peak_arr_ion[i,2] = en_peak 
        peak_arr_ion[i,3] = ef_peak 
        peak_arr_ion[i,4] = tau_peak
        peak_arr_ion[i,5] = min_tau_p

        peak_arr_no_ion[i,0] = beam_en_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr_no_ion[i])
        peak_arr_no_ion[i, 1] = p_peak 
        peak_arr_no_ion[i,2] = en_peak 
        peak_arr_no_ion[i,3] = ef_peak 
        peak_arr_no_ion[i,4] = tau_peak
        peak_arr_no_ion[i,5] = min_tau_p

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue", label="ion.")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,2]*1e9, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red", label="no ion.")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,2]*1e9, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue", label="ion.")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,3]*1e2, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red", label="no ion.")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,3]*1e2, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue", label="ion.")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,1], color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red", label="no ion.")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,1], color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()  

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated UV pulse durations ", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Mininum pulse duration (fs)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,4]*1e15, color="blue", label="ion.")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,4]*1e15, color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,4]*1e15, color="red", label="no ion.")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,4]*1e15, color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"min_tau_vs_beam_power.png"),dpi=1000)
    if show: plt.show()   

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("Gas: "+gas+"; CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated minimum pulse duration pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Mininum pulse duration pressure (bar)")
    ax.scatter(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,5], color="blue", label="ion.")
    ax.plot(peak_arr_ion[:,0]*1e6, peak_arr_ion[:,5], color="blue")
    ax.scatter(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,5], color="red", label="no ion.")
    ax.plot(peak_arr_no_ion[:,0]*1e6, peak_arr_no_ion[:,5], color="red")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()   

# plot UV energy and THG conversion efficiency for different beam power
def plot_phi_scan(sup_dir, gas, beam_en, ion, kerr, dens_mod):

    # create output directory
    if ion=="true":
        ion_string ="ion"
    else:    
        ion_string="no-ion"    

    out_path = os.path.join(out_dir, "phi_scan_"+gas+"_{0}mW_".format(beam_en*1e6)+kerr+"_"+ion_string+"_"+dens_mod)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # get data 
    path_arr, phi_arr, I = phi_comp(sup_dir, gas, beam_en,ion, kerr, dens_mod)

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = phi_arr / np.max(phi_arr) 

    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _, _, _ = get_data(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr, en_arr*1e9, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr, _, _ = get_data(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr, ef_arr*1e2, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: pulse duration vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Simulated UV pulse durations", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, tau_arr, _ = get_data(path_arr[i])
        plt.scatter(p_arr, tau_arr*1e15, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr, tau_arr*1e15, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

    # PLOT 4: z_peak vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, _, z_peak_arr = get_data(path_arr[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=cmap(cidx[i]), label="{0:.3f} rad".format(phi_arr[i]))
        plt.plot(p_arr,z_peak_arr*1e3, color=cmap(cidx[i]))

    plt.legend()
    plt.savefig(os.path.join(out_path,"z_peak.png"),dpi=1000)
    if show: plt.show()

    # PLOTS 5-9: peak investigation 
    peak_arr = np.empty((len(path_arr),6))

    for i in np.arange(len(path_arr)):
        peak_arr[i,0] = phi_arr[i]
        p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr[i])
        peak_arr[i, 1] = p_peak 
        peak_arr[i,2] = en_peak 
        peak_arr[i,3] = ef_peak
        peak_arr[i,4] = tau_peak
        peak_arr[i,5] = min_tau_p 

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    plt.xlabel("CEP (rad)")
    plt.ylabel("Peak UV energy (nJ)")
    plt.scatter(peak_arr[:,0], peak_arr[:,2]*1e9, color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,2]*1e9, color="blue")

    plt.savefig(os.path.join(out_path,"peak_en_vs_phi.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    plt.xlabel("CEP (rad)")
    plt.ylabel("Efficiency (%)")
    plt.scatter(peak_arr[:,0], peak_arr[:,3]*1e2, color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,3]*1e2, color="blue")

    plt.savefig(os.path.join(out_path,"peak_ef_vs_phi.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    plt.xlabel("CEP (rad)")
    plt.ylabel("Saturation pressure (bar)")
    plt.scatter(peak_arr[:,0], peak_arr[:,1], color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,1], color="blue")

    plt.savefig(os.path.join(out_path,"peak_p_vs_phi.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated minimum UV pulse durations", fontsize=16)
    plt.xlabel("CEP (rad)")
    plt.ylabel("Minimum pulse duration (fs)")
    plt.scatter(peak_arr[:,0], peak_arr[:,4]*1e15, color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,4]*1e15, color="blue")

    plt.savefig(os.path.join(out_path,"min_tau_vs_phi.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=[7.04, 5.28])
    plt.title("Gas: "+gas+"; beam power: {0}mW (I= {1:.1f}PW/cm2); ".format(beam_en*1e6, I*1e-15)+" response function: "+kerr+"; ionisation: "+ion, fontsize=10)
    plt.suptitle("Simulated minimum pulse duration pressures", fontsize=16)
    plt.xlabel("CEP (rad)")
    plt.ylabel("Minimum pulse duration pressure (bar)")
    plt.scatter(peak_arr[:,0], peak_arr[:,5], color="blue")
    plt.plot(peak_arr[:,0], peak_arr[:,5], color="blue")

    plt.savefig(os.path.join(out_path,"min_tau_p_vs_phi.png"),dpi=1000)
    if show: plt.show()

# single power, single dens model gas comparison 
def plot_gas_comp_singleP(sup_dir,beam_en,dens_mod):

    # get data 
    path_arr, gas_arr = gas_comp_singleP(sup_dir,beam_en,dens_mod)

    # get trivia
    _, phi, _, ion, kerr, _, _, I = get_params(path_arr[0])
   
    # create output directory
    if ion=="true":
        ion_string ="ion"
    else:    
        ion_string="no-ion"    

    out_path = os.path.join(out_dir, "gas_comp_"+str(beam_en*1e3)+"mW_"+str(phi)+"rad_"+kerr+"_"+ion_string+"_"+dens_mod)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

   # PLOT 1: energy vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Beam power: {1}mW (I= {2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _,_,_ = get_data(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, en_arr*1e9, color=colour_cycle[i])

    plt.legend()
    plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Beam power: {1}mW (I= {2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr, _, _ = get_data(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, ef_arr*1e2, color=colour_cycle[i])

    plt.legend()
    plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: pulse duration vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Simulated UV pulse durations", fontsize=16)
    plt.title("Beam power: {1}mW (I= {2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, tau_arr,_ = get_data(path_arr[i])
        plt.scatter(p_arr, tau_arr*1e15, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, tau_arr*1e15, color=colour_cycle[i])

    plt.legend()
    plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

    # PLOT 4: z_peak vs pressure 
    plt.figure(figsize=[7.04, 5.28]) 
    plt.subplots_adjust(top=0.85)
    plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title("Beam power: {1}mW (I= {2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, _,z_peak_arr = get_data(path_arr[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, z_peak_arr*1e3, color=colour_cycle[i])

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_positions.png"),dpi=1000)
    if show: plt.show()

# multi power, single dens model gas comparison 
def plot_gas_comp_multiP(sup_dir,dens_mod, excluded_gases):
    
    # get data 
    data, w0, tau = gas_comp_multiP(sup_dir,dens_mod,excluded_gases)
    gas_arr = data[:,0]

    # get trivia
    _, phi, _, ion, kerr, _, _, I = get_params(data[0,1][0])
   
    # create output directory
    if ion=="true":
        ion_string ="ion"
    else:    
        ion_string="no-ion"    

    out_path = os.path.join(out_dir, "gas_comp_multiP_"+str(phi)+"rad_"+kerr+"_"+ion_string+"_"+dens_mod)
    if not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # auxiliary functions (beam_en in mW, I in PW/cm^2):
    def p2i(beam_en):
        return beam_en *1e-6 / (np.pi * (w0*1e2)**2 * tau ) * 1e-15

    def i2p(I):
        return I * (np.pi * (w0*1e2)**2 * tau ) *1e15 *1e6    
    
    # peak investigations 
    peak_data = np.empty((len(gas_arr), 2), dtype="object")

    for j in np.arange(len(gas_arr)):

        path_arr = data[j, 1]
        beam_en_arr = data[j, 2]
        peak_arr = np.empty((len(path_arr),6))

        for i in np.arange(len(path_arr)):
            peak_arr[i,0] = beam_en_arr[i]
            p_peak, en_peak, ef_peak, tau_peak, min_tau_p = get_peaks(path_arr[i])
            peak_arr[i, 1] = p_peak 
            peak_arr[i,2] = en_peak 
            peak_arr[i,3] = ef_peak 
            peak_arr[i,4] = tau_peak 
            peak_arr[i,5] = min_tau_p

        peak_data[j, 0] = gas_arr[j]
        peak_data[j, 1] = peak_arr    

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color=colour_cycle[k])

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")

    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color=colour_cycle[k])
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color=colour_cycle[k])
    
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated UV pulse durations", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration (fs)")
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color=colour_cycle[k])
    
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    plt.savefig(os.path.join(out_path,"tau_min_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=[7.04, 5.28])
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.suptitle("Simulated minimal pulse duration pressure", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration pressure (bar)")
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,5], color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,5], color=colour_cycle[k])
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

# ---------- EXEC --------------------------------------
if single:
    plot_single(single_dir, n)
else:
   #for gas in ["Ar", "Ne", "He", "N2", "Kr", "N2O", "Xe"]: 
        #plot_beamP_scan(sup_dir, gas, 0.0, "true", "f", "coms") 
        #plot_beamP_scan(sup_dir, gas, 0.0, "true", "f", "grad") 
        #plot_multiP_singleg_model_comp(sup_dir, gas)
   
    placeholder = "_"
    plot_gas_comp_singleP(sup_dir,400e-3,"coms")
    #plot_gas_comp_multiP(sup_dir,"coms", [])
    #plot_gas_comp_multiP(sup_dir,"grad", [])
    #plot_singlePg_model_comp(sup_dir, "N2O", 400*1e-3)
