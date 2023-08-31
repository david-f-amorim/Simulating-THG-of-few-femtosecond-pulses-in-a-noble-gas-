import numpy as np
import matplotlib.pyplot as plt 
import os 

"""
THE PURPOSE OF THIS FILE IS TO PRODUCE RELEVANT SPECTRAL PLOTS FOR THE 
MANUSCRIPT. FILE NAMES WILL HAVE TO BE CHANGED MANUALLY BETWEEN PLOTS 
AND UNUSED CODE CHUNKS COMMENTED OUT.
"""

# ---------- INPUT/OUTPUT HANDLING  --------------------------------------
in_dir =  "manuscript_spectra"   # directory from which to read input files 
out_dir = "manuscript_spectra"     # directory to store output files 

# ---------- FORMAT PLOTS  --------------------------------------
disable_latex = False # toggle LaTeX rendering
show_title = False     # toggle showing title 
norm = True           # toggle normalisation
use_pdf = True        # if true: use pdf and png; else: use png only

if disable_latex == False : plt.rcParams["text.usetex"] = True   # enable LaTeX rendering
plt.rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
plt.rcParams["font.family"] = "STIXGeneral" # use LateX font for text
plt.rcParams["font.size"] = 16 # set standard font size 
fig_dim = [2 * 3.14961,2* 2.3622075] # for 8cm width ; double for 16cm width

# ---------- SELECT PLOTS -------------------------------

plot1 = False  # exp. and sim.
plot2 = True # old v. new 
plot3 = False # chirp 
plot4 = False # CEP
plot5 = False # 10 best 
# ---------- PLOTS  --------------------------------------

if plot1:
    # * * * PLOT 1: SHOW MEASURED SPECTRUM [WITH OR WITHOUT OVERLAYED SIMULATION]
    #                for manuscript: Ne 400mW 2.0bar; Ar 150mW 0.4bar

    file_measured = "Ar_150mW_0.4bar.txt"     # file name of measured UV spectrum 
    file_sim      =  "spectrum_Ar.txt"    # file name of simulated UV spectrum
    file_out      = "spec_comp_Ar_150mW_2.5scale_2.0bar"     # name of output file (no file ending!)

    overlay = True         # if True: overlay simulated spectrum 

    spec_measured = np.loadtxt(os.path.join(in_dir, file_measured))   # read in measured spectrum 
    if overlay: spec_sim = np.loadtxt(os.path.join(in_dir, file_sim)) # read in simulated spectrum 

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV spectra")
    plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
    plt.xlabel("$\lambda$ (nm)") 

    if norm:
        if overlay==False:
            spec_measured[:,1]=spec_measured[:,1]/np.max(spec_measured[:,1])
        else:
            spec_measured[:,1]=spec_measured[:,1]/ np.max(spec_measured)
            spec_sim[:,1]=spec_sim[:,1]/ np.max(spec_sim)    

    plt.plot(spec_measured[:,0]*1e9, spec_measured[:,1], color="red", label="exp.")
    plt.xlim(200, 360)
    if overlay:
        plt.plot(spec_sim[:,0]*1e9, spec_sim[:,1], color="blue", label="sim.")
        plt.legend()

    if use_pdf:
         plt.savefig(os.path.join(out_dir,file_out+".pdf")) 
    plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.show()

if plot2:
    # * * * PLOT 2: SHOW SIMULATIONS FOR OLD VERSUS NEW CELL
    #                for manuscript: Ne 400mW 2.0bar; Ar 150mW 0.4bar

    file_old = "spectrum_Ar_old.txt"     # file name of simulated UV spectrum (old cell)
    file_new = "spectrum_Ar_new.txt"     # file name of simulated UV spectrum (new cell)
    file_grad ="spectrum_Ar.txt" 
    file_out = "Ar_old_v_new_v_grad"     # name of output file (no file ending!)

    spec_new = np.loadtxt(os.path.join(in_dir, file_new))   # read in "new" spectrum
    spec_old = np.loadtxt(os.path.join(in_dir, file_old))   # read in "old" spectrum
    spec_grad = np.loadtxt(os.path.join(in_dir, file_grad))

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV spectra")
    plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
    plt.xlabel("$\lambda$ (nm)") 

    if norm:
        max_val = max(np.max(spec_new), np.max(spec_old),np.max(spec_grad))
        spec_new[:,1]=spec_new[:,1]/ max_val
        spec_old[:,1]=spec_old[:,1]/ max_val
        spec_grad[:,1]=spec_grad[:,1]/ max_val   

    plt.plot(spec_old[:,0]*1e9, spec_old[:,1], color="blue", label="old cell")
    plt.plot(spec_new[:,0]*1e9, spec_new[:,1], color="red", label="new chip")
    plt.plot(spec_grad[:,0]*1e9, spec_grad[:,1], color="black", label="gradient")

    plt.xlim(200,360)
    plt.legend()

    if use_pdf:
         plt.savefig(os.path.join(out_dir,file_out+".pdf")) 
         plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.show()

if plot3:
    # * * * PLOT 3: SHOW SIMULATIONS FOR DIFFERENT GVD VALUES
    #                for manuscript: Ne 400mW 2.0bar; Ar 150mW 0.4bar;
    #                                GVD: 0, +11fs^2, -11fs^2

    file_no_chirp = "spectrum_Ar.txt"     # file name of simulated UV spectrum (no chirp)
    file_pos_chirp="spectrum_Ar_pos_chirp.txt"      # file name of simulated UV spectrum (pos. chirp)
    file_neg_chirp="spectrum_Ar_neg_chirp.txt"      # file name of simulated UV spectrum (neg. chirp)
    file_out = "Ar_chirp"          # name of output file (no file ending!)

    spec_no_chirp = np.loadtxt(os.path.join(in_dir, file_no_chirp))   
    spec_pos_chirp = np.loadtxt(os.path.join(in_dir, file_pos_chirp))   
    spec_neg_chirp = np.loadtxt(os.path.join(in_dir, file_neg_chirp))   

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV spectra")
    plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
    plt.xlabel("$\lambda$ (nm)") 
    plt.xlim(200,360)

    if norm:
            max_val = max(np.max(spec_no_chirp), np.max(spec_pos_chirp), np.max(spec_neg_chirp))
            spec_no_chirp[:,1]=spec_no_chirp[:,1]/ max_val 
            spec_pos_chirp[:,1]=spec_pos_chirp[:,1]/ max_val
            spec_neg_chirp[:,1]=spec_neg_chirp[:,1]/ max_val   

    plt.plot(spec_pos_chirp[:,0]*1e9, spec_pos_chirp[:,1], color="red", label="+11.0fs$^2$")
    plt.plot(spec_neg_chirp[:,0]*1e9, spec_neg_chirp[:,1], color="blue", label="$-$11.0fs$^2$")
    plt.plot(spec_no_chirp[:,0]*1e9, spec_no_chirp[:,1], color="black", label="0.0fs$^2$")
    plt.legend()

    if use_pdf:
         plt.savefig(os.path.join(out_dir,file_out+".pdf")) 
         plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.show()

if plot4:
    # * * * PLOT 4: SHOW SIMULATIONS FOR DIFFERENT CEP VALUES
    #                for manuscript: Ne 400mW 2.0bar; Ar 150mW 0.4bar;
    #                                CEP: 0, pi/4, pi/2, 3pi/4, pi

    file_CEP_0 = ""     # file name of simulated UV spectrum (CEP 0)
    file_CEP_pi_4 = ""     # file name of simulated UV spectrum (CEP pi/4)
    file_CEP_pi_2 = ""     # file name of simulated UV spectrum (CEP pi/2)
    file_CEP_3pi_4 = ""     # file name of simulated UV spectrum (CEP 3pi/4)
    file_CEP_pi = ""     # file name of simulated UV spectrum (CEP pi)

    file_out = ""          # name of output file (no file ending!)

    spec_CEP_0 = np.loadtxt(os.path.join(in_dir, file_CEP_0))   
    spec_CEP_pi_4 = np.loadtxt(os.path.join(in_dir, file_CEP_pi_4))
    spec_CEP_pi_2 = np.loadtxt(os.path.join(in_dir, file_CEP_pi_2))
    spec_CEP_3pi_4 = np.loadtxt(os.path.join(in_dir, file_CEP_3pi_4))
    spec_CEP_pi = np.loadtxt(os.path.join(in_dir, file_CEP_pi))

    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV spectra")
    plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
    plt.xlabel("$\lambda$ (nm)") 

    if norm:
            max_val = max(np.max(spec_CEP_0), np.max(spec_CEP_pi_4), np.max(spec_CEP_pi_2), np.max(spec_CEP_3pi_4), np.max(spec_CEP_pi))
            spec_CEP_0[:,1]=spec_CEP_0[:,1]/ max_val
            spec_CEP_pi_4[:,1]=spec_CEP_pi_4[:,1]/ max_val
            spec_CEP_pi_2[:,1]=spec_CEP_pi_2[:,1]/ max_val 
            spec_CEP_3pi_4[:,1]=spec_CEP_3pi_4[:,1]/ max_val
            spec_CEP_pi[:,1]=spec_CEP_pi[:,1]/ max_val 

    cmap = plt.get_cmap("viridis")
    cval = np.linspace(0,1,5)

    plt.plot(spec_CEP_0[:,0]*1e9, spec_CEP_0[:,1], color=cmap(cval[0]), label="0")
    plt.plot(spec_CEP_pi_4[:,0]*1e9, spec_CEP_pi_4[:,1], color=cmap(cval[1]), label="$\pi$/4")
    plt.plot(spec_CEP_pi_2[:,0]*1e9, spec_CEP_pi_2[:,1], color=cmap(cval[2]), label="$\pi$/2")
    plt.plot(spec_CEP_3pi_4[:,0]*1e9, spec_CEP_3pi_4[:,1], color=cmap(cval[3]), label="3$\pi$/4")
    plt.plot(spec_CEP_pi[:,0]*1e9, spec_CEP_pi[:,1], color=cmap(cval[4]), label="$\pi$")
    plt.legend()

    if use_pdf:
         plt.savefig(os.path.join(out_dir,file_out+".pdf")) 
    plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.show()

if plot5:
    # * * * PLOT 5: SHOW THE NINE BEST SPECTRA
    #        

    file1 = ""
    file2 = ""
    file3 = ""
    file4 = ""
    file5 = ""
    file6 = ""
    file7 = ""
    file8 = ""
    file9 = ""

    file_out = ""          # name of output file (no file ending!)

    spec1 = np.loadtxt(os.path.join(in_dir, file1))  
    spec2 = np.loadtxt(os.path.join(in_dir, file2))  
    spec3 = np.loadtxt(os.path.join(in_dir, file3))  
    spec4 = np.loadtxt(os.path.join(in_dir, file4))  
    spec5 = np.loadtxt(os.path.join(in_dir, file5))  
    spec6 = np.loadtxt(os.path.join(in_dir, file6))  
    spec7 = np.loadtxt(os.path.join(in_dir, file7))  
    spec8 = np.loadtxt(os.path.join(in_dir, file8))   
    spec9 = np.loadtxt(os.path.join(in_dir, file9))  
   
    plt.figure(figsize=fig_dim) 
    plt.subplots_adjust(top=0.9, bottom=0.14)
    if show_title: plt.title("Simulated UV spectra")
    plt.ylabel("I (arb. units)" if norm==False else "I (norm.)")
    plt.xlabel("$\lambda$ (nm)") 

    if norm:
            max_val = max(np.max(spec1), np.max(spec2), np.max(spec3), np.max(spec4), np.max(spec5),np.max(spec6), np.max(spec7), np.max(spec8), np.max(spec9))
            spec1[:,1]=spec1[:,1]/ max_val
            spec2[:,1]=spec2[:,1]/ max_val
            spec3[:,1]=spec3[:,1]/ max_val
            spec4[:,1]=spec4[:,1]/ max_val
            spec5[:,1]=spec5[:,1]/ max_val
            spec6[:,1]=spec6[:,1]/ max_val
            spec7[:,1]=spec7[:,1]/ max_val
            spec8[:,1]=spec8[:,1]/ max_val
            spec9[:,1]=spec9[:,1]/ max_val

    colour_cycle = ["blue", "grey", "black", "red", "purple", "green", "cyan", "orange", "deeppink"]

    plt.plot(spec1[:,0]*1e9, spec1[:,1], color=colour_cycle[0], label="a")
    plt.plot(spec2[:,0]*1e9, spec2[:,1], color=colour_cycle[1], label="b")
    plt.plot(spec3[:,0]*1e9, spec3[:,1], color=colour_cycle[2], label="c")
    plt.plot(spec4[:,0]*1e9, spec4[:,1], color=colour_cycle[3], label="d")
    plt.plot(spec5[:,0]*1e9, spec5[:,1], color=colour_cycle[4], label="e")
    plt.plot(spec6[:,0]*1e9, spec6[:,1], color=colour_cycle[5], label="f")
    plt.plot(spec7[:,0]*1e9, spec7[:,1], color=colour_cycle[6], label="g")
    plt.plot(spec8[:,0]*1e9, spec8[:,1], color=colour_cycle[7], label="h")
    plt.plot(spec9[:,0]*1e9, spec9[:,1], color=colour_cycle[8], label="i")
    
    plt.legend(ncol=3)

    if use_pdf:
         plt.savefig(os.path.join(out_dir,file_out+".pdf")) 
    plt.savefig(os.path.join(out_dir,file_out+".png"),dpi=1000)
    plt.show()