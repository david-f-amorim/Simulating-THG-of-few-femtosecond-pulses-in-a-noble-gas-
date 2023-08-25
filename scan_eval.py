import numpy as np
import matplotlib.pyplot as plt
import os  
import json

# ---------- QUICK SETTINGS --------------------------------------

single = False  # if True: process output from a single pressure scan; if False: process output from a set of different pressure scans ("multi-scan")

n    = 8   # maximum number of overlayed spectra in one plot (to avoid clutter)
show = True # if True: open plots in matplotlib GUI; if False: don't show plots (directly write to file) 
save = True # if True: saves plots 

old = False # if True: skip some features that were not available with some older scan data; only use if error occurs when old == False

show_title = True
norm = True 
comp_exp = True
disable_latex = False 

shift_sim = "no" # either: "no", "offset", or "factor"
set_shift = None # if 

# ---------- MULTI-SCAN SETTINGS --------------------------------------



# ---------- INPUT/OUTPUT HANDLING --------------------------------------

single_dir = "parameter_scans\\gas_scans\\scan_200.0mW_Ne_0.0rad_f_ion_coms" # path to pressure scan directory to process if single=True (Note: output will be written to same directory)
exp_file = "raw_input\\Ne_200mW_IR.txt"

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
                    elif (key in params_dict) and (value !=params_dict.get(key)):
                        cond *=0 
                    else:
                        raise SystemExit("Error: argument '{0}' not recognised; kwargs should be: 'gas', 'dens_mod', 'tau', 'lam0', 'w0', 'CEP', 'IR_energy', 'ion', 'propz', 'GVD', 'thickness', 'ion_mod', 'IR_int'".format(key))    

            if cond: 
                IR_energy = float(params_arr[6])
                IR_int    = float(params_arr[12])

                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                IR_energy_arr = np.append(IR_energy_arr, IR_energy)
                IR_int_arr = np.append(IR_int_arr, IR_int)

    if len(path_arr)==0:
        raise SystemExit("Error: No files were matched to the specified selection.")            
    
    path_arr = path_arr[IR_energy_arr.argsort()] 
    IR_int_arr = IR_int_arr[IR_energy_arr.argsort()]
    IR_energy_arr = IR_energy_arr[IR_energy_arr.argsort()]
             

    return path_arr, IR_energy_arr, IR_int_arr


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

# get paths to two pressure scan files that differ only in 
# one parameter (set by second_var):

def two_comparison_get(sup_dir, second_var, **kwargs):

    path_arr =np.array([])
    kwargs_list = ['gas', 'dens_mod', 'tau', 'lam0', 'w0', 'CEP', 'IR_energy', 'ion', 'propz', 'GVD', 'thickness', 'ion_mod', 'IR_int']
    val_arr = np.array([])

    if second_var[0] in kwargs:
        raise SystemExit("Error: '{0}' cannot be used as 'second_var' and be in 'kwargs' simultaneously".format(second_var[0]))
    if (second_var[0] in kwargs_list)== False:
        raise SystemExit("Error: argument '{0}' not recognised; kwargs should be: 'gas', 'dens_mod', 'tau', 'lam0', 'w0', 'CEP', 'IR_energy', 'ion', 'propz', 'GVD', 'thickness', 'ion_mod', 'IR_int'".format(second_var[0]))    

    for dir in os.listdir(sup_dir):

        if os.path.isdir(os.path.join(sup_dir,dir)):
            
            _, params_dict = get_params(os.path.join(sup_dir,dir))

            cond = 1

            for key, value in kwargs.items():
                    if (key in params_dict)==False:
                        raise SystemExit("Error: argument '{0}' not recognised; kwargs should be: 'gas', 'dens_mod', 'tau', 'lam0', 'w0', 'CEP', 'IR_energy', 'ion', 'propz', 'GVD', 'thickness', 'ion_mod', 'IR_int'".format(key))    
                    if (key in params_dict) and (value !=params_dict.get(key)):
                        cond *=0    
                        

            if  (second_var[1]==params_dict.get(second_var[0])):
                cond *=1
                val = 1
            elif (second_var[2]==params_dict.get(second_var[0])):
                cond*=1
                val = 2
            else:
                cond *=0

            if cond: 
                path_arr = np.append(path_arr, os.path.join(sup_dir,dir))
                val_arr = np.append(val_arr, val)

    if len(path_arr)==0:
        raise SystemExit("Error: No files were matched to the specified selection.")
    elif len(path_arr)==1:
        raise SystemExit("Error: Only one file was matched to the specified selection.")
    elif len(path_arr) > 2:
        raise SystemExit("Error: Too many files were matched to the specified selection.")

    path_arr = path_arr[val_arr.argsort()]

    return path_arr

# get paths to 


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

        data[i, 0] = gas_arr[i]
        data[i, 1] = path_arr 
        data[i, 2] = beam_en_arr
    
    return data, w0, tau     

# ---------- FUNCTIONS THAT EXTRACT DATA FROM A SINGLE PRESSURE SCAN DIRECTORY --------------------------------------

# extract parameter values of the pressure scan directory at 'path' and returns array and dict of parameters
def get_params(path):

    params = np.loadtxt(os.path.join(path,"params.txt"),dtype="str", delimiter="=", comments="#")

    gas  = params[0,1][1:]
    dens_mod =params[2,1][1:]
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

    if  old==False:
        arr    = np.loadtxt(os.path.join(path,"energy_efficiency_time_zpeak.txt"))
    else:
        arr    = np.loadtxt(os.path.join(path,"energy_efficiency.txt"))
    arr    = arr[arr[:, 0].argsort()]
    p_arr  = arr[:,0] 
    UVen_arr = arr[:,1]
    ef_arr = arr[:,2]
    if old==False: 
        tau_arr= arr[:,3]
        zpeak_arr= arr[:,4]
    else:
        tau_arr= None
        zpeak_arr= None    

    return p_arr, UVen_arr, ef_arr, tau_arr, zpeak_arr 

# extract peak energy, peak efficiency, minimum pulse duration, pressure at peak energy, and pressure at minimum pulse duration of the pressure scan directory at 'path';
# returns array 
def get_peaks(path):

    p_arr, en_arr, ef_arr, tau_arr, _ = get_data(path)

    peak_en = np.max(en_arr)
    peak_ef = np.max(ef_arr)

    if old:
        min_tau = None
        min_tau_pres = None
    else:
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

# plot UV energy, THG conversion efficiency, pulse duration, UV peak position, and UV spectra for a single scan 
def plot_single(single_dir, n=n):

    # set plot formatting 
    if disable_latex == False : plt.rcParams["text.usetex"] = True   # enable LaTeX renadering
    plt.rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
    plt.rcParams["font.family"] = "STIXGeneral" # use LateX font for text
    plt.rcParams["font.size"] = 16 # set standard font size 
    fig_dim = [2 * 3.14961,2* 2.3622075] # for 8cm width ; double for 16cm width

    # get data 
    p_arr, en_arr, ef_arr, tau_arr, zpeak_arr = get_data(single_dir)
    peak_arr = get_peaks(single_dir)
    
    # import measured energy data (should be: first col pressure in bar; other cols energy in nJ)
    if comp_exp:
        measured_arr = np.loadtxt(exp_file, skiprows=0, delimiter=";", dtype="float")
        p_measured = measured_arr[:, 0]
        en_measured = np.mean(measured_arr[:,1:], axis=1)
        
        if shift_sim:
            if shift_sim == "offset" :

                if set_shift == None:
                    dif = peak_arr[3] - p_measured[np.where(en_measured == np.max(en_measured) )][0]
                else: dif = set_shift

                p_arr -= dif
            elif shift_sim == "factor" :
                if set_shift == None:
                    fac = peak_arr[3] / p_measured[np.where(en_measured == np.max(en_measured) )][0]
                else: fac = set_shift
                p_arr /= fac
            
    # plots
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV energies")
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, en_arr*1e9, color="blue")
    plt.scatter(p_arr, en_arr*1e9, color="blue", label= ("Peak: {0:.1f}nJ at {1}bar".format(peak_arr[0]*1e9, peak_arr[3]) if comp_exp == False else "sim."))
    if comp_exp:
        plt.plot(p_measured, en_measured, color="red")
        plt.scatter(p_measured, en_measured, color="red", label="exp.")
        plt.xlim(0, max(p_measured)+0.1)

    plt.legend(loc="upper right")

    if save: plt.savefig(os.path.join(single_dir,"energies.png"),dpi=1000)
    if show: plt.show()

    if comp_exp==False:
    
        plt.figure(figsize=fig_dim) 
        if show_title: plt.title("Simulated THG efficiencies")
        plt.subplots_adjust(top=0.9, bottom=0.14, left=0.16)
        plt.ylabel('Efficiency (\%)')
        plt.xlabel("Central pressure (bar)")
        plt.plot(p_arr, ef_arr*1e2, color="blue")
        plt.scatter(p_arr, ef_arr*1e2, color="blue", label="Peak: {0:.2f}\% at {1}bar".format(peak_arr[1]*1e2, peak_arr[3]))
        plt.legend()

        if save: plt.savefig(os.path.join(single_dir,"efficiencies.png"),dpi=1000)
        if show: plt.show()

        if old==False:
            plt.figure(figsize=fig_dim) 
            plt.subplots_adjust(top=0.9, bottom=0.14)
            if show_title: plt.title("Simulated pulse durations")
            plt.ylabel("Pulse duration (fs)")
            plt.xlabel("Central pressure (bar)")
            plt.plot(p_arr, tau_arr*1e15, color="blue")
            plt.scatter(p_arr, tau_arr*1e15, color="blue", label="Minimum: {0:.2f}fs at {1}bar".format(peak_arr[2]*1e15, peak_arr[4]))
            plt.legend()

            if save: plt.savefig(os.path.join(single_dir,"durations.png"),dpi=1000)
            if show: plt.show()

            plt.figure(figsize=fig_dim) 
            plt.subplots_adjust(top=0.9, bottom=0.14)
            if show_title: plt.title("Position of peak UV energy")
            plt.ylabel("Position (mm)")
            plt.xlabel("Central pressure (bar)")
            plt.plot(p_arr, zpeak_arr*1e3, color="blue")
            plt.scatter(p_arr, zpeak_arr*1e3, color="blue", label="Maximum: {0:.2f}mm at {1}bar".format(np.max(zpeak_arr)*1e3,p_arr[np.where(zpeak_arr == np.max(zpeak_arr) )][0]))
            plt.legend()

            if save: plt.savefig(os.path.join(single_dir,"z_peak.png"),dpi=1000)
            if show: plt.show()

        plt.figure(figsize=fig_dim) 
        plt.subplots_adjust(top=0.9, bottom=0.14)
        if show_title: plt.title("Simulated UV spectra")
        plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
        plt.xlabel("$\lambda$ (nm)")

        N = len(p_arr)
        cmap = plt.get_cmap("viridis")
        data, p_cut = get_spectra(single_dir, n)
        cidx = p_cut / np.max(p_cut) 
        
        if N<=n:
            n=N

        if norm: 
            maxI = 0
            for i in np.arange(n):
                if np.max(data[i,2]) > maxI:
                    maxI = np.max(data[i,2])    

        for i in np.arange(n):
            plt.plot(data[i,1]*1e9, data[i,2] if norm==False else data[i,2]/maxI, color=cmap(cidx[i]), label="{0}bar".format(data[i,0]))
        plt.legend(loc="upper right")

        if save: plt.savefig(os.path.join(single_dir,"spectra.png"),dpi=1000)
        if show: plt.show()

# plot UV energy, THG conversion efficiency, pulse duration, UV peak position, for two scans, differing in second_var 
def plot_double(scan_dir, second_var,**kwargs):

    # get paths to two files
    paths = two_comparison_get(scan_dir, second_var,**kwargs)

    # get data 
    p_arr, UVen_arr, ef_arr, tau_arr, zpeak_arr = get_data(paths[0])
    peak_arr = get_peaks(paths[0])

    p_arr2, UVen_arr2, ef_arr2, tau_arr2, zpeak_arr2 = get_data(paths[1])
    peak_arr2 = get_peaks(paths[1])

    # set labels
    if second_var[0]=="gas": 
        label1 = second_var[1]
        label2 = second_var[2]
    elif second_var[0]=="dens_mod":
        if second_var[1]=="coms":
            label1="COMSOL"
            label2="grad."
        else:
            label1="grad."
            label2="COMSOL"
    elif second_var[0]=="ion":
        if second_var[1]=="true":
            label1="ion."
            label2="no ion."
        else:
            label1="no ion."
            label2="ion." 
    elif second_var[0]=="ion_mod":
        if second_var[1]=="ADK":
            label1="ADK"
            label2="PPT"
        else:
            label1="PPT"
            label2="ADK"
    else:
        label1=second_var[1]
        label2=second_var[2]                                   

    # set output string
    out_dir_str=""
    for key, val in kwargs.items():
        out_dir_str += key+"="+str(val)+"_"
    out_path = os.path.join(out_dir, "two_comparison_"+out_dir_str) 
    if save and not os.path.isdir(out_path): 
        os.mkdir(out_path)    

    # set plot formatting 
    if disable_latex == False : plt.rcParams["text.usetex"] = True   # enable LaTeX renadering
    plt.rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
    plt.rcParams["font.family"] = "STIXGeneral" # use LateX font for text
    plt.rcParams["font.size"] = 16 # set standard font size 
    fig_dim = [2 * 3.14961,2* 2.3622075] # for 8cm width ; double for 16cm width

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.15)
    if show_title: plt.title("Simulated UV energies")
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, UVen_arr*1e9, color="blue")
    plt.scatter(p_arr, UVen_arr*1e9, color="blue", label=label1)
    plt.plot(p_arr2, UVen_arr2*1e9, color="red")
    plt.scatter(p_arr2, UVen_arr2*1e9, color="red", label=label2)
    plt.legend(loc="upper left")

    if save: plt.savefig(os.path.join(out_path,"energies.png"),dpi=1000)
    if show: plt.show()

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.15)
    if show_title: plt.title("Simulated THG efficiencies")
    plt.ylabel("Efficiency (\%)")
    plt.xlabel("Central pressure (bar)")
    plt.plot(p_arr, ef_arr*1e2, color="blue")
    plt.scatter(p_arr, ef_arr*1e2, color="blue", label=label1)
    plt.plot(p_arr2, ef_arr2*1e2, color="red")
    plt.scatter(p_arr2, ef_arr2*1e2, color="red", label=label2)
    plt.legend(loc="upper left")

    if save: plt.savefig(os.path.join(out_path,"efficiencies.png"),dpi=1000)
    if show: plt.show()

    if old==False:
        plt.figure(figsize=fig_dim) 
        plt.subplots_adjust(top=0.9, bottom=0.15)
        if show_title: plt.title("Simulated pulse durations")
        plt.ylabel("Pulse duration (fs)")
        plt.xlabel("Central pressure (bar)")
        plt.plot(p_arr, tau_arr*1e15, color="blue")
        plt.scatter(p_arr, tau_arr*1e15, color="blue", label=label1)
        plt.plot(p_arr2, tau_arr2*1e15, color="red")
        plt.scatter(p_arr2, tau_arr2*1e15, color="red", label=label2)
        plt.legend()

        if save: plt.savefig(os.path.join(out_path,"durations.png"),dpi=1000)
        if show: plt.show()

        plt.figure(figsize=fig_dim) 
        plt.subplots_adjust(top=0.9, bottom=0.15)
        if show_title: plt.title("Position of peak UV energy")
        plt.ylabel("Position (mm)")
        plt.xlabel("Central pressure (bar)")
        plt.plot(p_arr, zpeak_arr*1e3, color="blue")
        plt.scatter(p_arr, zpeak_arr*1e3, color="blue", label=label1)
        plt.plot(p_arr2, zpeak_arr2*1e3, color="red")
        plt.scatter(p_arr2, zpeak_arr2*1e3, color="red", label=label2)
        plt.legend()

        if save: plt.savefig(os.path.join(out_path,"z_peak.png"),dpi=1000)
        if show: plt.show()

# plot UV energy and THG conversion efficiency for different beam power
# NOTE: a second_var can be overlayed (meant mainly for ion_mod, dens_mod,ion, and gas)
#       for second_var: both options should have the same number of relevant files!!
def plot_multipower(sup_dir, second_var=(None, None, None), **kwargs):

    # check if second_var is valid (if given)
    if (second_var[0] != None) and (second_var[1] != None) and (second_var[0] != None):
        if second_var[0] in kwargs:
           raise SystemExit("Error: '{0}' cannot be used as 'second_var' and be in 'kwargs' simultaneously".format(second_var[0]))
        else:
            bool_sec_var = True
    else:
        bool_sec_var = False   

    # set output directory name 
    out_dir_str=""
    for key, val in kwargs.items():
        out_dir_str += key+"="+str(val)+"_"     

    # get data and params
    if bool_sec_var:
        kwargs2 = kwargs.copy()
        kwargs[second_var[0]]=second_var[1]
        kwargs2[second_var[0]]=second_var[2]
        path_arr, IR_energy_arr, IR_int_arr = power_comparison_get(sup_dir,**kwargs)
        path_arr2, IR_energy_arr2, IR_int_arr2 = power_comparison_get(sup_dir,**kwargs2)

        if len(path_arr) != len(path_arr2):
            raise SystemExit("Error: the selected files with '{0}={1}' cannot be matched up with the files with '{0}={2}'. The selected files should only differ in '{0}' and there should be a corresponding '{0}={1}' file for each '{0}={2}' file.".format(second_var[0],second_var[1],second_var[2]))

        # set labels
        if second_var[0]=="gas": 
            label1 = second_var[1]
            label2 = second_var[2]
        elif second_var[0]=="dens_mod":
            if second_var[1]=="coms":
                label1="COMSOL"
                label2="grad."
            else:
                label1="grad."
                label2="COMSOL"
        elif second_var[0]=="ion":
            if second_var[1]=="true":
                label1="ion."
                label2="no ion."
            else:
                label1="no ion."
                label2="ion." 
        elif second_var[0]=="ion_mod":
            if second_var[1]=="ADK":
                label1="ADK"
                label2="PPT"
            else:
                label1="PPT"
                label2="ADK"
        else:
            label1=second_var[1]
            label2=second_var[2]
    else:
        path_arr, IR_energy_arr, IR_int_arr = power_comparison_get(sup_dir,**kwargs)    

    N = len(path_arr)
    
    p_arr, UVen_arr, ef_arr, tau_arr, zpeak_arr = np.empty(N, dtype="object"),np.empty(N, dtype="object"), np.empty(N, dtype="object"), np.empty(N, dtype="object"),np.empty(N, dtype="object")
    
    if bool_sec_var:
        p_arr2, UVen_arr2, ef_arr2, tau_arr2, zpeak_arr2 = np.empty(N, dtype="object"),np.empty(N, dtype="object"), np.empty(N, dtype="object"), np.empty(N, dtype="object"),np.empty(N, dtype="object")
    
    for i in np.arange(N):
       p_arr[i], UVen_arr[i], ef_arr[i], tau_arr[i], zpeak_arr[i] = get_data(path_arr[i])
       if bool_sec_var:
            p_arr2[i], UVen_arr2[i], ef_arr2[i], tau_arr2[i], zpeak_arr2[i] = get_data(path_arr2[i])

    params, params_dict = get_params(path_arr[0])
    w0, tau = float(params[4]), float(params[2])

    # set output directory 
    out_path = os.path.join(out_dir, "beam_power_comparison_"+out_dir_str) 
    if save and not os.path.isdir(out_path): 
        os.mkdir(out_path)

    # auxiliary functions (beam_power in mW, I in PW/cm^2):
    def p2i(beam_power):
        return beam_power *1e-6 / (np.pi * (w0*1e2)**2 * tau ) * 1e-15

    def i2p(I):
        return I * (np.pi * (w0*1e2)**2 * tau ) *1e15 *1e6
    
    # set plot formatting 
    if disable_latex == False : plt.rcParams["text.usetex"] = True   # enable LaTeX renadering
    plt.rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
    plt.rcParams["font.family"] = "STIXGeneral" # use LateX font for text
    plt.rcParams["font.size"] = 16 # set standard font size 
    fig_dim = [2 * 3.14961,2* 2.3622075] # for 8cm width ; double for 16cm width

    # PLOT 1: energy vs pressure 
    cmap = plt.get_cmap("viridis")
    cidx = IR_energy_arr / np.max(IR_energy_arr) 

    plt.figure(figsize=fig_dim) 
    if show_title: plt.title("Simulated UV energies")
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(N):

        if bool_sec_var:
            plt.scatter(p_arr[i], UVen_arr[i]*1e9, color=cmap(cidx[i]), label=label1)
            plt.plot(p_arr[i], UVen_arr[i]*1e9, color=cmap(cidx[i]))
            plt.scatter(p_arr2[i], UVen_arr2[i]*1e9, color=cmap(cidx[i]),marker="+", label=label2)
            plt.plot(p_arr2[i], UVen_arr2[i]*1e9, ls="--", color=cmap(cidx[i]))
        else:
            plt.scatter(p_arr[i], UVen_arr[i]*1e9, color=cmap(cidx[i]), label="{0}mW ({1:.1f}PW/cm^2)".format(IR_energy_arr[i]*1e6, 1e-15*IR_int_arr[i]))
            plt.plot(p_arr[i], UVen_arr[i]*1e9, color=cmap(cidx[i]))

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=fig_dim) 
    if show_title: plt.title("Simulated THG efficiencies")
    plt.ylabel("Efficiency (\%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        
        if bool_sec_var:
            plt.scatter(p_arr[i], ef_arr[i]*1e2, color=cmap(cidx[i]), label=label1)
            plt.plot(p_arr[i], ef_arr[i]*1e2, color=cmap(cidx[i]))
            plt.scatter(p_arr2[i],ef_arr2[i]*1e2, color=cmap(cidx[i]),marker="+", label=label2)
            plt.plot(p_arr2[i], ef_arr2[i]*1e2, ls="--", color=cmap(cidx[i]))
        else:
            plt.scatter(p_arr[i], ef_arr[i]*1e2, color=cmap(cidx[i]), label="{0}mW ({1:.1f}PW/cm^2)".format(IR_energy_arr[i]*1e6, 1e-15*IR_int_arr[i]))
            plt.plot(p_arr[i], ef_arr[i]*1e2, color=cmap(cidx[i]))

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: pulse duration vs pressure
    if old==False:  
        plt.figure(figsize=fig_dim) 
        if show_title: plt.title("Simulated UV pulse durations")
        plt.ylabel("Pulse duration (fs)")
        plt.xlabel("Central pressure (bar)")
        
        for i in np.arange(len(path_arr)):
            if bool_sec_var:
                plt.scatter(p_arr[i], tau_arr[i]*1e15, color=cmap(cidx[i]), label=label1)
                plt.plot(p_arr[i], tau_arr[i]*1e15, color=cmap(cidx[i]))
                plt.scatter(p_arr2[i],tau_arr2[i]*1e15,marker="+", color=cmap(cidx[i]), label=label2)
                plt.plot(p_arr2[i], tau_arr2[i]*1e15,ls="--", color=cmap(cidx[i]))
            else:
                plt.scatter(p_arr[i], tau_arr[i]*1e15, color=cmap(cidx[i]), label="{0}mW ({1:.1f}PW/cm^2)".format(IR_energy_arr[i]*1e6, 1e-15*IR_int_arr[i]))
                plt.plot(p_arr[i], tau_arr[i]*1e15, color=cmap(cidx[i]))

        plt.legend()
        if save: plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
        if show: plt.show()

    # PLOT 4: z_peak vs pressure 
    if old==False:
        plt.figure(figsize=fig_dim) 
        if show_title: plt.title("Position of peak UV energy")
        plt.ylabel("Position (mm)")
        plt.xlabel("Central pressure (bar)")
        
        for i in np.arange(len(path_arr)):

            if bool_sec_var:
                plt.scatter(p_arr[i], zpeak_arr[i]*1e3, color=cmap(cidx[i]), label=label1)
                plt.plot(p_arr[i], zpeak_arr[i]*1e3, color=cmap(cidx[i]))
                plt.scatter(p_arr2[i],zpeak_arr2[i]*1e3,marker="+", color=cmap(cidx[i]), label=label2)
                plt.plot(p_arr2[i], zpeak_arr2[i]*1e3,ls="--", color=cmap(cidx[i]))
            else:
                plt.scatter(p_arr[i], zpeak_arr[i]*1e3, color=cmap(cidx[i]), label="{0}mW ({1:.1f}PW/cm^2)".format(IR_energy_arr[i]*1e6, 1e-15*IR_int_arr[i]))
                plt.plot(p_arr[i], zpeak_arr[i]*1e3, color=cmap(cidx[i]))

        plt.legend()
        if save: plt.savefig(os.path.join(out_path,"peak_positions.png"),dpi=1000)
        if show: plt.show()

    # PLOTS 4-7: peak investigation 
    peak_arr = np.empty((len(path_arr),6))
    if bool_sec_var:
        peak_arr2 = np.empty((len(path_arr),6))

    for i in np.arange(len(path_arr)):
        peak_arr[i,0] = IR_energy_arr[i]
        peak_en, peak_ef, min_tau, peak_en_pres, min_tau_pres = get_peaks(path_arr[i])
        peak_arr[i, 1] = peak_en_pres 
        peak_arr[i,2] = peak_en 
        peak_arr[i,3] = peak_ef 
        peak_arr[i,4] = min_tau 
        peak_arr[i,5] =  min_tau_pres

        if bool_sec_var:
            peak_arr2[i,0] = IR_energy_arr2[i]
            peak_en, peak_ef, min_tau, peak_en_pres, min_tau_pres = get_peaks(path_arr2[i])
            peak_arr2[i,1] = peak_en_pres 
            peak_arr2[i,2] = peak_en 
            peak_arr2[i,3] = peak_ef 
            peak_arr2[i,4] = min_tau 
            peak_arr2[i,5] =  min_tau_pres

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    if show_title: plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    if bool_sec_var:
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue", label="{0}={1}".format(second_var[0], second_var[1]))
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")
        ax.scatter(peak_arr2[:,0]*1e6, peak_arr2[:,2]*1e9, color="red", label="{0}={1}".format(second_var[0], second_var[2]))
        ax.plot(peak_arr2[:,0]*1e6, peak_arr2[:,2]*1e9, color="red")
        plt.legend()
    else:    
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color="blue")

    if save: plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    if show_title: plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")

    if bool_sec_var:
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue", label="{0}={1}".format(second_var[0], second_var[1]))
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
        ax.scatter(peak_arr2[:,0]*1e6, peak_arr2[:,3]*1e2, color="red", label="{0}={1}".format(second_var[0], second_var[2]))
        ax.plot(peak_arr2[:,0]*1e6, peak_arr2[:,3]*1e2, color="red")
        plt.legend()
    else:
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color="blue")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    if save: plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title(title_str, fontsize=10)
    if show_title: plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    if bool_sec_var:
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue", label="{0}={1}".format(second_var[0], second_var[1]))
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
        ax.scatter(peak_arr2[:,0]*1e6, peak_arr2[:,1], color="red", label="{0}={1}".format(second_var[0], second_var[2]))
        ax.plot(peak_arr2[:,0]*1e6, peak_arr2[:,1], color="red")
        plt.legend()
    else:
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color="blue")
    
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    if save: plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    if old==False:
        fig, ax = plt.subplots(figsize=fig_dim)
        plt.subplots_adjust(top=0.8, right=0.9)
        ax.set_title(title_str, fontsize=10)
        if show_title: plt.suptitle("Simulated UV pulse durations", fontsize=16)
        ax.set_xlabel("Beam power (mW)")
        ax.set_ylabel("Minimum pulse duration (fs)")
        if bool_sec_var:
            ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue", label="{0}={1}".format(second_var[0], second_var[1]))
            ax.plot(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue")
            ax.scatter(peak_arr2[:,0]*1e6, peak_arr2[:,4]*1e15, color="red", label="{0}={1}".format(second_var[0], second_var[2]))
            ax.plot(peak_arr2[:,0]*1e6, peak_arr2[:,4]*1e15, color="red")
            plt.legend()
        else:
            ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue")
            ax.plot(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color="blue")
        
        secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
        secax.set_xlabel('Peak intensity (PW/cm^2)')

        if save: plt.savefig(os.path.join(out_path,"tau_min_vs_beam_power.png"),dpi=1000)
        if show: plt.show()

    if old==False:
        fig, ax = plt.subplots(figsize=fig_dim)
        plt.subplots_adjust(top=0.8, right=0.9)
        ax.set_title(title_str, fontsize=10)
        if show_title: plt.suptitle("Simulated minimal pulse duration pressure", fontsize=16)
        ax.set_xlabel("Beam power (mW)")
        ax.set_ylabel("Minimum pulse duration pressure (bar)")
        if bool_sec_var:
            ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue", label="{0}={1}".format(second_var[0], second_var[1]))
            ax.plot(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue")
            ax.scatter(peak_arr2[:,0]*1e6, peak_arr2[:,5], color="red", label="{0}={1}".format(second_var[0], second_var[2]))
            ax.plot(peak_arr2[:,0]*1e6, peak_arr2[:,5], color="red")
            plt.legend()
        else:
            ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue")
            ax.plot(peak_arr[:,0]*1e6, peak_arr[:,5], color="blue")
        
        secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
        secax.set_xlabel('Peak intensity (PW/cm^2)')

        if save: plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
        if show: plt.show()

        # save kwargs and sup_dir to file
        params_dict["sup_dir"] = sup_dir

        if save: 
            with open(os.path.join(out_path,"kwargs.txt"), "w") as file:
                json.dump(params_dict, file)

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
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.85)
    if show_title: plt.suptitle("Simulated UV energies", fontsize=16)
    plt.title("Beam power: {1}mW ({2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Energy (nJ)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, en_arr, _,_,_ = get_data(path_arr[i])
        plt.scatter(p_arr, en_arr*1e9, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, en_arr*1e9, color=colour_cycle[i])

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"UV_energies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 2: efficiency vs pressure 
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.85)
    if show_title: plt.suptitle("Simulated THG efficiencies", fontsize=16)
    plt.title("Beam power: {1}mW ({2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Efficiency (%)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, ef_arr, _, _ = get_data(path_arr[i])
        plt.scatter(p_arr, ef_arr*1e2, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, ef_arr*1e2, color=colour_cycle[i])

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"THG_efficiencies.png"),dpi=1000)
    if show: plt.show()

    # PLOT 3: pulse duration vs pressure 
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.85)
    if show_title: plt.suptitle("Simulated UV pulse durations", fontsize=16)
    plt.title("Beam power: {1}mW ({2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Pulse duration (fs)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, tau_arr,_ = get_data(path_arr[i])
        plt.scatter(p_arr, tau_arr*1e15, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, tau_arr*1e15, color=colour_cycle[i])

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"pulse_durations.png"),dpi=1000)
    if show: plt.show()

    # PLOT 4: z_peak vs pressure 
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.85)
    if show_title: plt.suptitle("Position of peak UV energy", fontsize=16)
    plt.title("Beam power: {1}mW ({2:.1f}PW/cm^2); CEP: {0:.2f}rad; ".format(phi,beam_en*1e3, I*1e-15 )+"\n response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    plt.ylabel("Position (mm)")
    plt.xlabel("Central pressure (bar)")
    
    for i in np.arange(len(path_arr)):
        p_arr, _, _, _,z_peak_arr = get_data(path_arr[i])
        plt.scatter(p_arr, z_peak_arr*1e3, color=colour_cycle[i], label=gas_arr[i])
        plt.plot(p_arr, z_peak_arr*1e3, color=colour_cycle[i])

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"peak_positions.png"),dpi=1000)
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

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    if show_title: plt.suptitle("Simulated peak UV energies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Peak UV energy (nJ)")
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,2]*1e9, color=colour_cycle[k])

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"peak_en_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    if show_title: plt.suptitle("Simulated peak UV efficiencies", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Efficiency (%)")

    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,3]*1e2, color=colour_cycle[k])
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"peak_ef_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    if show_title: plt.suptitle("Simulated saturation pressures", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Saturation pressure (bar)")
    
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,1], color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,1], color=colour_cycle[k])
    
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')

    plt.legend()
    if save: plt.savefig(os.path.join(out_path,"peak_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    if show_title: plt.suptitle("Simulated UV pulse durations", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration (fs)")
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,4]*1e15, color=colour_cycle[k])
    
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    if save: plt.savefig(os.path.join(out_path,"tau_min_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

    fig, ax = plt.subplots(figsize=fig_dim)
    plt.subplots_adjust(top=0.8, right=0.9)
    ax.set_title("CEP: {0:.2f}rad; ".format(phi)+"response function: "+kerr+"; ionisation: "+ion+"; model: "+dens_mod, fontsize=10)
    if show_title: plt.suptitle("Simulated minimal pulse duration pressure", fontsize=16)
    ax.set_xlabel("Beam power (mW)")
    ax.set_ylabel("Minimum pulse duration pressure (bar)")
    for k in np.arange(len(gas_arr)):
        peak_arr = peak_data[k, 1]
        ax.scatter(peak_arr[:,0]*1e6, peak_arr[:,5], color=colour_cycle[k], label=peak_data[k,0])
        ax.plot(peak_arr[:,0]*1e6, peak_arr[:,5], color=colour_cycle[k])
    secax = ax.secondary_xaxis('top', functions=(p2i, i2p))
    secax.set_xlabel('Peak intensity (PW/cm^2)')
    plt.legend()

    if save: plt.savefig(os.path.join(out_path,"min_tau_p_vs_beam_power.png"),dpi=1000)
    if show: plt.show()

# ---------- EXEC --------------------------------------
if single:
    plot_single(single_dir, n)
else:  
    plot_double(sup_dir, ("dens_mod", "coms", "grad"),IR_energy=75e-6, gas="Ar")
