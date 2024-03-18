using  Luna
import Luna.PhysData: wlfreq    
import FFTW                    
import Luna: Hankel  
import NumericalIntegration: integrate, SimpsonEven          
import Dates                   
using  DelimitedFiles
using  LaTeXStrings 
       
# ----------------- QUICK SETTINGS -------------------------------
p_scan = false               # if true, pressure scan is executed, with pressure range set by the variable "pres_arr" below; if false, a single run is simulated
save = true                  # if true, save output plots, run parameters and UV spectrum
show = false                 # if true, opens plots in GUI after run 
txt_only = false             # if true, no plots are produced (only UV spectrum and simulation parameters are written to file)

read_IR = false               # if true: read input IR pulse [time domain] from file; if false: use Gaussian approximation 
read_ρ  = false               # if true: read gas density profile from file; if false: use pressure gradient approximation 
fluo_profile = false          # if true: use gas density profile from fluorescence profile; only necessary if read_ρ = true
use_gaussfit = false         # if true: use gas density profile from gaussian fit to fluorescence profile; only necessary if read_ρ = true and fluo_profile = true
read_UV = false              # if true: read in measured UV spectrum and overlay on simulated results   
only_r0 = false              # if false: E-field is integrated over r for time domain analysis, intensity analysis, ...

IR_spec = false              # if true: read measured input IR spectrum from file and overlay 
show_IR = false              # if true and "read_IR" is true: overlay measured time-domain input pulse on plots (barely ever relevant; used to check that transform to frequency domain works)
IR_spec_exp = false          # if true: read input IR spectrometer spectrum from file and overlay (barely ever relevant; used to check that FROG spectrum is accurate)

save_UV_temp = true          # if true: write time-domain output UV pulse to file (has no effect for pressure scans)
save_UV_phase = true
save_phase = true
save_UV_en = true	    # if true: write z-dependant UV energy in txt output file

# ------------------ SET PHYSICAL PARAMETERS ------------------------

gas = :Ne           # gas type in cell 
pres = 1e-3 #5.0          # central gas pressure [bar] (for single run; not relevant for pressure scans)
p_ed = 1e-3         # edge gas pressure [bar] 
p_const = false     # if true: set constant pressure profile P==(pres,pres,pres) ; if false: set simple gradient: P==(p_ed, pres, p_ed); (only relevant when read_ρ==false)
τ = 3e-15 #5e-15           # FWHM pulse duration of IR pulse [s] (only relevant when read_IR==false)
λ0 = 260e-9 #800e-9         # central wavelength of IR pulse [m]
w0 = 200e-6 #65e-6          # beam waist of IR pulse [m]
CEP = 0.0           # carrier-envelope phase of IR pulse (ϕ0) [rad] (only relevant when read_IR==true)                                 
IRenergy = 1000e-9 #400e-6    # IR pulse energy [J] (related to beam power via 1kHz repetition rate)                                                           
#L = 5.58e-3            # propagation distance (cell length) [m]
L = 25.0e-3            # propagation distance (cell length) [m]

propz = -L/2        # propagation distance from the waist [m], i.e. beam focus position  (NOTE: always specified in a coordinate system where the cell starts at 0 and ends at L!)
z_vals =L .* [1.0/4, 1.0/3, 3.0/4, 1.0]     # points along the cell at which to investigate beam evolution [m] (NOTE: always specified in a coordinate system where the cell starts at 0 and ends at L!)

λ_lims = (100e-9, 1000e-9)       # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (100e-9, 360e-9)     # wavelength limits of UV region of interest [m,m]
λ_rangeIR = (600e-9, 1000e-9)    # wavelength limits of IR region of interest [m,m]

material = :SiO2     # material of the optics the IR beam propagates through before THG (for chirp compensation)
thickness= 0         # thickness of the material [m] (for chirp compensation)

ϕs = [0,0,0*1e-30,0]     # Taylor-series coefficients of initial spectral phase  [s^n] (used to introduce additional chirp)

ion = false         # if true: enable ionisation response, if false: disable ionisation 
ion_model="PPT"    # set to "ADK" or "PPT" (has no effect if ion==false); "ADK" is less accurate at low intensities but faster; "PPT" may crash at very high intensities

pres_arr = range(start= 2.0, stop= 6.0, step= 0.25)  # pressure range (only relevant for pressure scans)  [bar]


# ----------------- PLOT SETTINGS -------------------------------

show_title = false # if false: hides figure titles
norm = true        # if true: normalise intensities on figures
disable_latex = false # if true: disable latex rendering of plots (saves time but might result in some labels or title being displayed incorrectly)
use_pdf = false     # if true: save output as pdf; if false: use png

# ----------------- FILE HANDLING -------------------------------

# NOTE: see "./file_prepare.py" for relevant file content specifications

in_dir    = "input"                          # directory of input files

file_IR   = "IRpulse_old.dat"                    # name of IR input pulse file 
path_IR   = joinpath(in_dir, file_IR)        # sys. path to IR input pulse file 

file_UV   = "UVpulse.dat"                    # name of UV output pulse file 
path_UV   = joinpath(in_dir, file_UV)        # sys. path to UV output pulse file 

file_IR_spec = "IRspec.dat"                   # name of IR FROG input spectrum file 
path_IR_spec = joinpath(in_dir, file_IR_spec) # sys. path to IR FROG input spectrum file 

file_IR_spec_exp = "IRspec_exp.dat"                   # name of IR input spectrometer spectrum file 
path_IR_spec_exp = joinpath(in_dir, file_IR_spec_exp) # sys. path to IR input spectrometer spectrum file 

# ----------------- MAKE ARRANGEMENTS FOR PRESSURE SCAN -----------

if p_scan == true  
    txt_only = true        # ensures that no plots are generated as part of the scan (save time)
    save     = true        # ensures that output is saved 
end    

scan_dir = "scan_"*string(IRenergy*1e6)*"mW_"*string(gas)*"_"*string(round(CEP; digits=3))*"rad_"*string(ion ? "ion" : "no-ion")*"_"*string(read_ρ ? "coms" : "grad")  # name of scan output directory 

# ----------------- DEFINE MAIN FUNCTION -----------

function THG_main(pres=pres) 

    # ----------------- SET PRESSURE PROFILE ---------------------------

    if read_ρ == true 
        
        if fluo_profile == true
            
            if use_gaussfit == true
                file_ρ    = "dens_0.25bar_FLUOGAUSS.dat"    # name of density profile data file 
            else
                file_ρ    = "dens_0.25bar_FLUORESCENCE.dat"    # name of density profile data file
            end
            
            path_ρ    = joinpath(in_dir, file_ρ)                  # sys. path to density profile data file 
            
            max_ρ  = PhysData.density(gas, pres)#, T=293)         # returns the number density [m^-3] at pressure `pres` and temperature 293 K = 20 C.

            z_in = readdlm(path_ρ,' ', Float64, '\n')[:,1]        # read in spatial points 
            ρ_in = readdlm(path_ρ,' ', Float64, '\n')[:,2]        # read in fluorescence data 
            ρ_in = ρ_in/maximum(ρ_in)*max_ρ                       # normalize density profile to the maximum density defined by the variable "pres"

        else
            file_ρ    = "dens_$(pres)bar.dat"                     # name of density profile data file 
            path_ρ    = joinpath(in_dir, file_ρ)                  # sys. path to density profile data file 

            z_in = readdlm(path_ρ,' ', Float64, '\n')[:,1]        # read in spatial points 
            ρ_in = readdlm(path_ρ,' ', Float64, '\n')[:,2]        # read in density data 

        end
            
        dens = Maths.CSpline(z_in, ρ_in)                      # interpolate density function 
        γ = PhysData.sellmeier_gas(gas)                       # get polarisability from Sellmeier expansion

        coren(ω; z) = sqrt(1 + γ(wlfreq(ω)*1e6)*dens(z))      # calculate refractive index along the cell   

        L_total = maximum(z_in)                               # get total propagation distance         
        #print("\nmaximum(z_in)", maximum(z_in))
        #print("\nρ_in", ρ_in)
        
    else   
        Z= (0, L/2, L)                                        # define points for pressure gradient
        p_const ? P=(pres,pres,pres) : P= (p_ed, pres, p_ed)  # define values for pressure gradient (see definition "p_const" above)
        (coren,dens)=Capillary.gradient(gas,Z,P)              # gives n(ω; z) and ρ(z) from pressure profile   
    end  

    # ----------------- SET SIMULATION GRID ----------------------------

    R = 250e-6           # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  
    N = 1024              # sample size for Hankel tansform grid  (NOTE: this value was taken from an older version of the code; no justification! )           
    trange = 0.05e-12     # total extent of time window required [s] (NOTE: if this is too short, the range is extended automatically; this value was taken from an older version of the code; no justification!)                            

    if read_ρ==true 
        grid = Grid.RealGrid(maximum(z_in), λ0, λ_lims ,trange)   # set up time & space grid for read-in density case 
    else    
        grid = Grid.RealGrid(L, λ0, λ_lims ,trange)               # set up time & space grid for gradient approximation 
    end

    q = Hankel.QDHT(R, N, dim=2)                  # set up discrete Hankel transform matrix                        
    energyfun, _ = Fields.energyfuncs(grid, q)    # "energyfun" gives total energy in a field E(t)    

    # ----------------- SET NONLINEAR EFFECTS ----------------------------

    kerr = "f"                                                  # set nonlinear Kerr effect (NOTE: legacy from an older version! leave at "f")
    ionpot = PhysData.ionisation_potential(gas)                 # set gas ionisation potential   
    
    if ion_model=="ADK"
        ionrate = Ionisation.ionrate_fun!_ADK(gas)                  # set gas ionisation rate (ADK)
    elseif ion_model=="PPT"
        ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)        # set gas ionisation rate (PPT)
    end    
    
    linop = LinearOps.make_linop(grid, q, coren)                # generate linear operator for pulse-propagation equation
    normfun = NonlinearRHS.norm_radial(grid, q, coren)          # generate normalisation function for radial symmetry 

    # * * * SET KERR EFFECT AND/OR PLASMA FORMATION 
    if (kerr=="f") & (ion == true)                              # nonlinear response function with ionisation and Kerr effect
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                )  
    elseif (kerr=="f") & (ion == false)                        # nonlinear response function without ionisation, just Kerr effect
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                )    
    end    

    # ----------------- SET INPUT FIELD ----------------------------                                       

    if read_IR == true 

        t_in = readdlm(path_IR,' ', Float64, '\n')[:,1]                   # read in timing data
        I_in = readdlm(path_IR,' ', Float64, '\n')[:,2]                   # read in intensity data (time domain)

        beam_spline    = Maths.CSpline(t_in, I_in)                        # interpolate temporal beam intensity envelope  
        beam_spline_fs = Maths.CSpline(t_in*1e15, I_in)                   # interpolate temporal beam intensity envelope with t in fs

        function beam_profile(t,r::AbstractVector)
        beam_spline.(t) .* Maths.gauss.(r, w0/2)'                         # model spatial beam profile as Gaussian  
        end

        if read_ρ==false 
            inputs = Fields.SpatioTemporalField(λ0, IRenergy, CEP, 0.0,           # model input beam based off measured data (coordinate system with cell between 0 and L)
            (t, r) -> beam_profile(t, r), propz)  
        else 
            inputs = Fields.SpatioTemporalField(λ0, IRenergy, CEP, 0.0,           # model input beam based off measured data (shifted coordinate system due to regions outside cell)
            (t, r) -> beam_profile(t, r), propz+L_total/2)
        end 

    else 
        if read_ρ==false 
            inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ,                   # model input beam as Gaussian 
                    energy=IRenergy, w0=w0, propz=propz)
        else 
            inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ,                   # model input beam as Gaussian (shifted coordinate system due to regions outside cell)
            energy=IRenergy, w0=w0, propz=propz+L_total/2)
        end                
    end            

    # ---------------- SET UP SIMULATION -------------------

    Eω, transform, FT = Luna.setup(grid, q, dens, normfun, responses, inputs)  # set up propagation 

    # ----------------- ADD CHIRP ----------------------------

    Fields.prop_material!(Eω, grid, material, thickness, λ0)              # account for chirp due to mirrors, etc. 
    Fields.prop_taylor!(Eω, grid, ϕs, λ0)                                 # add chirp using spectral phase Taylor coefficients ϕs

    # ----------------- RUN SIMULATION ----------------------------

    output = Output.MemoryOutput(0, grid.zmax, 201)                        # configure output (NOTE: these values were taken from a previous version of the code)
    Luna.run(Eω, grid, linop, transform, FT, output)                       # run simulation  

    # ----------------- PROCESS RESULTS ----------------------------

    # * * * EXTRACT GRID PARAMS:  
    ω = grid.ω                   # sampled angular frequency values [rad/s]
    f = ω / 2π                   # sampled linear frequency values [Hz]
    λ = ω .\ (2π * PhysData.c)   # sampled wavelengths [m]         
    t = grid.t                   # sampled points in time [s]                                          
    zout = output.data["z"]      # sampled points along z [m] 

    # * * * DEFINE USEFUL GRID INDICES:
    ωhighUVidx    = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeUV[1]))         # index X such that ω[X]=ω_maxUV 
    ωlowUVidx     = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeUV[2]))         # index X such that ω[X]=ω_minUV
    ωhighIRidx    = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeIR[1]))         # index X such that ω[X]=ω_maxIR 
    ωlowIRidx     = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeIR[2]))         # index X such that ω[X]=ω_minIR

    λhighidx    = argmin(abs.(λ .- λ_rangeUV[1]))                         # index X such that λ[X]=λ_maxUV 
    λlowidx     = argmin(abs.(λ .-λ_rangeUV[2]))                          # index X such that λ[X]=λ_minUV

    # * * * EXTRACT FIELD AMPLITUDES:
    Eout = output.data["Eω"]                         # Hankel-transformed amplitude in frequency domain  
    Erout = (q \ Eout)                               # Real-space amplitude in frequency domain at r ≠ 0   (NOTE:"\" represents inverse Hankel transform) 
    
    if only_r0
	Er0 = dropdims(Hankel.onaxis(Eout, q), dims=2)
    else
        Er0 = zeros(ComplexF64, (size(Eout, 1), size(Eout, 3)))      # set up array for total (integrated along r) real-space amplitude in frequency domain 
    	for i = 1:size(Eout, 1), j = 1:size(Eout, 3)
            Er0[i,j] = Hankel.integrateK(Eout[i,:,j], q)             # integrate along r (technically: k) to obtain total real-space amplitude in frequency domain  
        end
    end

    Etout = FFTW.irfft(Erout, length(t), 1)     # time-domain real field amplitude at r≠0
    Et0 = FFTW.irfft(Er0, length(t),1)          # total time-domain real field amplitude across all radii
    Et = Maths.hilbert(Etout)                   # time-domain real field amplitude of envelope at r≠0

    # * * * FILTER FOR UV FIELD (r≠0):
    filter=FFTW.rfft(Etout, 1)        # set up filter array 
    filter[1:ωlowUVidx,:,:].=0;       # filters out ω < ω_min
    filter[ωhighUVidx:end,:,:].=0;    # filters out ω > ω_max 

    Etout_UV=FFTW.irfft(filter, length(t), 1)    # time-domain real field amplitude of UV pulse 
    Et_UV = Maths.hilbert(Etout_UV)              # time-domain real field amplitude of UV envelope

    # * * * FILTER FOR UV FIELD (across r)
    filter_all_r = FFTW.rfft(Et0, 1)          # set up filter array
    filter_all_r[1:ωlowUVidx,:].=0;           # filters out ω < ω_min
    filter_all_r[ωhighUVidx:end,:].=0;        # filters out ω > ω_max  

    Et0_UV =FFTW.irfft(filter_all_r, length(t), 1)     # total time-domain real field amplitude of UV pulse 
    It0_UV = abs2.(Et0_UV)                              # intensity of total UV pulse at 

    # * * * FILTER FOR IR FIELD (r≠0):
    filter_IR=FFTW.rfft(Etout, 1)        # set up filter array 
    filter_IR[1:ωlowIRidx,:,:].=0;       # filters out ω < ω_min
    filter_IR[ωhighIRidx:end,:,:].=0;    # filters out ω > ω_max 

    Etout_IR=FFTW.irfft(filter_IR, length(t), 1)    # time-domain real field amplitude of IR pulse
    Et_IR = Maths.hilbert(Etout_IR)                 # time-domain real field amplitude of IR envelope

    # * * * FILTER FOR IR FIELD (across r)
    filter_all_r_IR = FFTW.rfft(Et0, 1)          # set up filter array
    filter_all_r_IR[1:ωlowIRidx,:].=0;           # filters out ω < ω_min
    filter_all_r_IR[ωhighIRidx:end,:].=0;        # filters out ω > ω_max  

    Et0_IR =FFTW.irfft(filter_all_r_IR, length(t), 1)    # time-domain real field amplitude of IR pulse across r
    It0_IR = abs2.(Et0_IR)                               # intensity of total IR pulse across r

    # * * * GET FREQUENCY-RESOLVED UV TEMPORAL PROFILE AT OUTPUT 
    n = 15 # bin size 

    E_ωt_UV = zeros(ComplexF64,(length(ω[ωlowUVidx:n:ωhighUVidx]), length(t)))

    for i = 1:length(ω[ωlowUVidx:n:ωhighUVidx])
        Er0_filtered = Er0[:,end]
        ωi = range(ωlowUVidx,ωhighUVidx; step=n)[i]
        Er0_filtered[1:ωi] .=0 
        Er0_filtered[ωi+n:end] .= 0
        trans= Maths.hilbert(FFTW.irfft(Er0_filtered[:,end], length(t),1)) # for each UV frequency component at the output find temporal envelope
        for j = 1:length(t)
            E_ωt_UV[i,j] = trans[j]
        end 
    end  
    
    I_ωt_UV = abs2.(E_ωt_UV)

    # * * * EXTRACT INTENSITY ENVELOPES 
    It0_envelope = abs2.(Maths.hilbert(Et0))          # envelope modulating It0
    It0_UV_envelope = abs2.(Maths.hilbert(Et0_UV))    # envelope modulating It0_UV
    It0_IR_envelope = abs2.(Maths.hilbert(Et0_IR))    # envelope modulating It0_IR 

    # * * * CALCULATE PULSE DURATIONS
    τ_input = Maths.fwhm(t, It0_envelope[:,1])             # FWHM input pulse duration     [s]
    τ_UV    = Maths.fwhm(t, It0_UV_envelope[:,end])        # FWHM UV output pulse duration [s]

    # * * * CALCULATE PULSE ENERGIES:
    tot_pulse_en = zeros(length(zout))                   # set up array for total pulse energy 
    UV_pulse_en  = zeros(length(zout))                   # set up array for UV pulse energy 

    for i = 1:size(Etout, 3)
        tot_pulse_en[i] = energyfun(Etout[:, :, i])         # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
        UV_pulse_en[i]  = energyfun(Etout_UV[:, :, i])      # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
    end     

    η_THG = UV_pulse_en[end]/tot_pulse_en[1]            # THG efficiency: initial total pulse energy / final UV pulse energy 

    z_peak = zout[findmax(UV_pulse_en)[2]]             # z-coordinate of peak UV energy

    # * * * CONVERT TO INTENSITIES:
    Iωr = abs2.(Erout)                               #  intensity at r≠0 in frequency domain  [arbitrary units]
    Iω0 = abs2.(Er0)                                 #  intensity across all r in frequency domain  [arbitrary units]
    It0 = abs2.(Et0)                                 #  intensity across all r in time domain [arbitrary units]
    It = abs2.(Maths.hilbert(Etout))                 #  intensity of envelope at r≠0 in time domain [arbitrary units]

    # * * * FIND SPECTRAL PHASE  
    ϕω0 = angle.(Er0)                                # spectral phase across all r 
    ϕωr = angle.(Erout)                              # spectral phase at r≠0

    # * * * INTEGRATE FREQUENCY DOMAIN UV
    #       AND IR INTENSITIES 
    Iωr_UV = zeros((length(q.r),length(zout)))                   # set up arrays
    Iωr_IR  =zeros((length(q.r),length(zout)))
    Iω0_UV = zeros(length(zout))                   
    Iω0_IR  = zeros(length(zout))
    
    for i = 1:length(zout)

        for j=1:length(q.r)
            Iωr_UV[j,i] = abs2.(integrate(ω[ωlowUVidx:ωhighUVidx],Erout[ωlowUVidx:ωhighUVidx, j, i], SimpsonEven())) # frequency domain UV intensity (at r ≠0)
            Iωr_IR[j,i] = abs2.(integrate(ω[ωlowIRidx:ωhighIRidx],Erout[ωlowIRidx:ωhighIRidx, j, i], SimpsonEven()))  # frequency domain IR intensity (at r ≠ 0)
        end    

        Er0_int_UV = integrate(ω[ωlowUVidx:ωhighUVidx],Er0[ωlowUVidx:ωhighUVidx,  i] ,SimpsonEven()) # integrated UV field across all r
        Er0_int_IR = integrate(ω[ωlowIRidx:ωhighIRidx],Er0[ωlowIRidx:ωhighIRidx,  i] ,SimpsonEven()) # integrated IR field across all r

        Iω0_UV[i] = abs2.(Er0_int_UV) # frequency domain UV intensity (integrated along r)
        Iω0_IR[i] = abs2.(Er0_int_IR) # frequency domain IR intensity (integrated along r)
    end    
    
    # * * * PROCESS MEASURED DATA FROM FILES 
    if (IR_spec == true) & (txt_only == false)
        λ_IR_spec = readdlm(path_IR_spec,' ', Float64, '\n')[:,1]             # read in IR input wavelength data [FROG]
        I_IR_spec = readdlm(path_IR_spec,' ', Float64, '\n')[:,2]             # read in IR input spectral data [FROG]
    end

    if (IR_spec_exp == true) & (txt_only == false)
        λ_IR_spec_exp = readdlm(path_IR_spec_exp,' ', Float64, '\n')[:,1]     # read in IR input wavelength data [spectrometer]
        I_IR_spec_exp = readdlm(path_IR_spec_exp,' ', Float64, '\n')[:,2]     # read in IR input spectral data [spectrometer]
    end

    if (read_UV == true) & (txt_only == false) 
        λ_in    = readdlm(path_UV,' ', Float64, '\n')[:,1]                     # read in UV wavelength data
        I_in_UV = readdlm(path_UV,' ', Float64, '\n')[:,2]                     # read in UV intensity data 
    end    

    # ----------------- OUTPUT HANDLING -------------------------

    if isdir("output")== false   # create general directory to store output
        mkdir("output")
    end    

    if (save & !p_scan) == true              # create specific directory to store run output 
        out_path = joinpath("output","run_"*Dates.format(Dates.now(), "yyyy_mm_dd__HH_MM_SS"))
        mkdir(out_path)
    elseif p_scan == true  
        out_path = joinpath("output", scan_dir)
    end 

    # ----------------- PLOT RESULTS ----------------------------
    if txt_only == false 

        @eval import PyPlot: pygui, plt, PyDict, matplotlib 
        close("all")
        pygui(true)

        # set plot formatting 
        rcParams = PyDict(matplotlib."rcParams") # get rcParams 
        if disable_latex==false  rcParams["text.usetex"] = true end # enable LaTeX renadering
        rcParams["mathtext.fontset"] = "cm" # use LateX font for maths
        rcParams["font.family"] = "STIXGeneral" # use LateX font for text
        rcParams["font.size"] = 16 # set standard font size 
        fig_dim = 2* [3.14961, 2.3622075] # for 8cm width ; double for 16cm width  
        
        #+++++ PLOT 1:  IR and UV intensities as functions of z and r≠0
        plt.figure(figsize=fig_dim)
        if show_title plt.suptitle("Off-axis intensity of IR and UV beams") end
        plt.subplots_adjust(left=0.125, bottom=0.11, right=0.992, top=0.87, hspace=0.6)

        plt.subplot(2,1,1)
        plt.pcolormesh(zout*1e3, q.r*1e3, norm ? Maths.normbymax(Iωr_IR) : Iωr_IR)
        plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)" ))
        plt.ylabel("r (mm)")
        plt.xlabel("z (mm)")
        plt.title("IR beam")    

        plt.subplot(2,1,2)
        plt.pcolormesh(zout*1e3, q.r*1e3,norm ? Maths.normbymax(Iωr_UV) : Iωr_UV)
        plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)" ))
        plt.xlabel("z (mm)")
        plt.ylabel("r (mm)")
        plt.title("UV beam")

        if save==true
            if use_pdf == true
                plt.savefig(joinpath(out_path,"off-axis_intensity.pdf"))
            else 
                plt.savefig(joinpath(out_path,"off-axis_intensity.png"),dpi=1000)
            end
        end
        
        #+++++ PLOT 2:  IR and UV intensities as functions of z
        plt.figure(figsize=fig_dim)
        if show_title plt.suptitle("Total intensity of IR and UV beams") end
        plt.subplots_adjust(hspace=0.6)

        plt.subplot(2,1,1)
        plt.plot(zout*1e3,norm ? Maths.normbymax(Iω0_IR) : Iω0_IR, color="red")
        plt.title("IR beam")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")
        plt.xlabel("z (mm)")
        if norm==false plt.ticklabel_format(axis="y", style="scientific", scilimits=(0,0)) end
        
        plt.subplot(2,1,2)
        plt.plot(zout*1e3,  norm ? Maths.normbymax(Iω0_UV) : Iω0_UV, color="red")
        plt.title("UV beam")
        plt.xlabel("z (mm)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")
        if norm==false plt.ticklabel_format(axis="y", style="scientific", scilimits=(0,0)) end

        plt.tight_layout()

        if save==true
            if use_pdf == true
                plt.savefig(joinpath(out_path,"total_intensity.pdf"))
            else 
                plt.savefig(joinpath(out_path,"total_intensity.png"),dpi=1000)
            end
        end  
    
        #+++++ PLOT 3: gas number density and effective susceptibility along the cell 
        plt.figure(figsize=fig_dim)
        if show_title plt.title("Gas density profile") end
        
        plt.plot(zout*1e3, [dens(i) for i in zout] ./ PhysData.N_A, label=L"P_0="*"$(pres) bar", color="red")
        plt.ylabel(L"\rho"*" (mol/m"*L"^3"*")")
        plt.xlabel("z (mm)")
        plt.legend(loc="upper right")

        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"density.pdf"))
            else     
                plt.savefig(joinpath(out_path,"density.png"),dpi=1000)
            end 
        end 
        
        #+++++ PLOT 4:  linear  spectrum I(λ) at z=0 and z=L 
        plt.figure(figsize=fig_dim)
        if show_title plt.title("Full spectrum") end
        plt.plot(λ[2:end]*1e9,norm ? Maths.normbymax(Iω0[2:end,1]) : Iω0[2:end,1], label="z=$(round(zout[1]*1e3,digits=2))mm", color="grey")
        plt.plot(λ[2:end]*1e9,norm ? Iω0[2:end,end]/maximum(Iω0[2:end,1]) : Iω0[2:end,end], label="z=$(round(zout[end]*1e3,digits=2))mm", color="red")
        plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
        plt.xlabel(L"\lambda"*" (nm)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9,norm ? maximum(Iω0[λlowidx:λhighidx,end])/maximum(Iω0[2:end,1]) .*Maths.normbymax(I_in_UV) : maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9,norm ? Maths.normbymax(I_IR_spec) : maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec), color="green", ls="--", label="IR spectrum FROG") 
        end   
        
        if IR_spec_exp == true 
            plt.plot(λ_IR_spec_exp*1e9,norm ? Maths.normbymax(I_IR_spec_exp) : maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec_exp), color="blue", ls="--", label="IR spectrum spectrometer") 
        end 

        plt.legend()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"full_spectrum.pdf"))
            else     
                plt.savefig(joinpath(out_path,"full_spectrum.png"),dpi=1000)
            end 
        end 
        
        #+++++ PLOT 5:  UV only linear spectrum I(λ) at z=0 and z=L 
        plt.figure(figsize=fig_dim)
        if show_title plt.title("UV spectrum") end
        plt.plot(λ[λlowidx:λhighidx]*1e9,norm ? Iω0[λlowidx:λhighidx,1]/maximum(Iω0[λlowidx:λhighidx,end]) : Iω0[λlowidx:λhighidx,1], label="z=$(round(zout[1]*1e3,digits=2))mm", color="grey")
        plt.plot(λ[λlowidx:λhighidx]*1e9,norm ? Maths.normbymax(Iω0[λlowidx:λhighidx,end]) : Iω0[λlowidx:λhighidx,end], label="z=$(round(zout[end]*1e3,digits=2))mm", color="red")
        plt.xlim(λ_rangeUV[1]*1e9, λ_rangeUV[2]*1e9)
        plt.xlabel(L"\lambda"*" (nm)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9,norm ? Maths.normbymax(I_in_UV) : maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9,norm ? maximum(Iω0[λlowidx:λhighidx,1])/maximum(Iω0[λlowidx:λhighidx,end]) .*Maths.normbymax(I_IR_spec) : maximum(Iω0[λlowidx:λhighidx,1]).*Maths.normbymax(I_IR_spec), color="green", ls="--", label="IR spectrum FROG") 
        end

        if IR_spec_exp == true 
            plt.plot(λ_IR_spec_exp*1e9,norm ? maximum(Iω0[λlowidx:λhighidx,1]) /maximum(Iω0[λlowidx:λhighidx,end]) .*Maths.normbymax(I_IR_spec_exp) : maximum(Iω0[λlowidx:λhighidx,1]).*Maths.normbymax(I_IR_spec_exp), color="blue", ls="--",label="IR spectrum spectrometer") 
        end

        plt.legend()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"UV_spectrum.pdf"))
            else     
                plt.savefig(joinpath(out_path,"UV_spectrum.png"),dpi=1000)
            end 
        end 
        
	#+++++ PLOT 5a:  UV only linear spectrum I(λ) at z=0 
        plt.figure(figsize=fig_dim)
        if show_title plt.title("UV spectrum") end
        plt.plot(λ[λlowidx:λhighidx]*1e9,norm ? Maths.normbymax(Iω0[λlowidx:λhighidx,end]) : Iω0[λlowidx:λhighidx,end], label="z=$(round(zout[end]*1e3,digits=2))mm", color="purple")
        plt.fill_between(λ[λlowidx:λhighidx]*1e9, norm ? Maths.normbymax(Iω0[λlowidx:λhighidx,end]) : Iω0[λlowidx:λhighidx,end],0 * Maths.normbymax(Iω0[λlowidx:λhighidx,end]), color="purple", alpha=0.2)
        #plt.xlim(λ_rangeUV[1]*1e9, λ_rangeUV[2]*1e9)
        plt.xlim([180, 330])
        plt.ylim([0, 1.1])
        plt.xlabel(L"\lambda"*" (nm)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9,norm ? Maths.normbymax(I_in_UV) : maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9,norm ? maximum(Iω0[λlowidx:λhighidx,1])/maximum(Iω0[λlowidx:λhighidx,end]) .*Maths.normbymax(I_IR_spec) : maximum(Iω0[λlowidx:λhighidx,1]).*Maths.normbymax(I_IR_spec), color="green", ls="--") 
        end

        if IR_spec_exp == true 
            plt.plot(λ_IR_spec_exp*1e9,norm ? maximum(Iω0[λlowidx:λhighidx,1]) /maximum(Iω0[λlowidx:λhighidx,end]) .*Maths.normbymax(I_IR_spec_exp) : maximum(Iω0[λlowidx:λhighidx,1]).*Maths.normbymax(I_IR_spec_exp), color="blue", ls="--",label="IR spectrum spectrometer") 
        end

        plt.xticks([200, 250, 300])
        #plt.yticks([])
        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"UV_spectrum_nice.pdf"))
		plt.savefig(joinpath(out_path,"UV_spectrum_nice.png"),dpi=1000)
		plt.savefig(joinpath(out_path,"UV_spectrum_nice.eps"))
            else     
                plt.savefig(joinpath(out_path,"UV_spectrum_nice.png"),dpi=1000)
            end 
        end
        
        #+++++ PLOT 6:  log. spectrum I(λ) at z=0 and z=L 
        Iω0log = log10.(Iω0)

        plt.figure(figsize=fig_dim)
        if show_title plt.title("Logarithmic spectrum") end
        plt.plot(λ[2:end]*1e9,norm ? Maths.normbymax(Iω0log[2:end,1]) : Iω0log[2:end,1], label="z=$(round(zout[1]*1e3,digits=2))mm",  color="grey")
        plt.plot(λ[2:end]*1e9,norm ? Iω0log[2:end,end]/maximum(Iω0log[2:end,1]) : Iω0log[2:end,end], label="z=$(round(zout[end]*1e3,digits=2))mm", color="red")
        plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
        plt.xlabel(L"\lambda"*" (nm)")
        plt.ylabel(norm ? "log. I (norm.)" : "log. I (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9,norm ? log10.(abs.(maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV))) /maximum(Iω0log[2:end,1]) : log10.(abs.(maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV))), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9,norm ? log10.(abs.(maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec))) /maximum(Iω0log[2:end,1]) :  log10.(abs.(maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec))), color="green", ls="--", label="IR spectrum FROG") 
        end

        if IR_spec_exp == true 
            plt.plot(λ_IR_spec_exp*1e9,norm ?  log10.(abs.(maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec_exp)))  /maximum(Iω0log[2:end,1]) : log10.(abs.(maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec_exp))), color="blue", ls="--", label="IR spectrum spectrometer") 
        end

        plt.legend()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"log_spectrum.pdf"))
            else     
                plt.savefig(joinpath(out_path,"log_spectrum.png"),dpi=1000)
            end 
        end 

        #+++++ PLOT 7:  pulse energies and efficiency 
        plt.figure(figsize=fig_dim)
        if show_title plt.suptitle("Pulse energies ("*L"\eta_{THG}="*string(round(η_THG, digits=4) *100)*L"\% )") end
        plt.subplots_adjust(hspace=0.6, top=0.86)

        plt.subplot(2,1,1)
        plt.xlabel("z (mm)")
        plt.ylabel("E ("*L"\mu"*"J)")
        plt.plot(zout.*1e3, tot_pulse_en.*1e6, label=L"\Delta"*"E=-$(round(Int64, tot_pulse_en[1]*1e6-tot_pulse_en[end]*1e6))"*L"\mu"*"J", color="red")
        plt.title("Total pulse ")
        plt.legend()

        plt.subplot(2,1,2)
        plt.plot(zout.*1e3, UV_pulse_en.*1e9, label=L"\Delta"*"E=+$(round(Int64, UV_pulse_en[end]*1e9))nJ", color="red")
        plt.xlabel("z (mm)")
        plt.ylabel("E (nJ)")
        plt.title("UV pulse")
        plt.legend()

        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"pulse_energies.pdf"))
            else     
                plt.savefig(joinpath(out_path,"pulse_energies.png"),dpi=1000)
            end 
        end 

        if save_UV_en==true
            open(joinpath(out_path,"total_energy_evolution.txt"), "w") do file
                writedlm(file, zip(zout,tot_pulse_en))
            end
            open(joinpath(out_path,"UV_energy_evolution.txt"), "w") do file
                writedlm(file, zip(zout,UV_pulse_en))
            end
        end

        #+++++ PLOT 8: frequency evolution 
        plt.figure(figsize=fig_dim)
        if show_title plt.suptitle("Frequency evolution") end
        plt.pcolormesh(zout*1e3, f*1e-15,log10.(Maths.normbymax(Iω0)))   
        plt.clim(-6, 0)    
        plt.colorbar(label="log. I (arb. units)")
        plt.xlabel("z (mm)")
        plt.ylabel("f (PHz)")

        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"frequency_evolution.pdf"))
            else     
                plt.savefig(joinpath(out_path,"frequency_evolution.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 9: time-domain plot of input pulse 
        plt.figure(figsize=fig_dim) 
        if show_title plt.title("Time-domain representation of input pulse") end
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel(norm ? "I(z=0) (norm.)" : "I(z=0) (arb. units)")
        plt.plot(t*1e15,norm ? Maths.normbymax(It0[:,1]) : It0[:,1] , color="red", label=L"\tau_{FWHM}="*string(round(τ_input*1e15, digits=1) )*"fs")
        plt.plot(t*1e15,norm ? Maths.normbymax(It0_envelope[:,1]) : It0_envelope[:,1], color="black", ls="--")

        if (read_IR & show_IR )==true  # overlay measured input pulse 
            plt.plot(t_in*1e15,norm ? Maths.normbymax(I_in)  : maximum(It0_envelope[:,1]).*Maths.normbymax(I_in), color="green")
        end

        plt.legend(loc="upper right")
        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"time_domain_input.pdf"))
            else     
                plt.savefig(joinpath(out_path,"time_domain_input.png"),dpi=1000)
            end 
        end
        
        #+++++ PLOT 10: time-domain plot of UV output pulse 
        plt.figure(figsize=fig_dim) 
        if show_title plt.title("Time-domain representation of UV output pulse") end
        plt.xlabel("t (fs)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")
        plt.plot(t*1e15, norm ? Maths.normbymax(It0_UV[:,end]) : It0_UV[:,end] , color="red", label=L"\tau_{FWHM}="*string(round(τ_UV*1e15, digits=1) )*"fs")
        plt.plot(t*1e15,norm ? Maths.normbymax(It0_UV_envelope[:,end]) : It0_UV_envelope[:,end], color="black", linestyle="--")

        plt.legend(loc="upper right")
        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"time_domain_UV.pdf"))
		plt.savefig(joinpath(out_path,"time_domain_UV.png"),dpi=1000)
		plt.savefig(joinpath(out_path,"time_domain_UV.eps"))
            else     
                plt.savefig(joinpath(out_path,"time_domain_UV.png"),dpi=1000)
            end 
        end
        
        #+++++ PLOT 10a: time-domain plot of output IR pulse 
         
        plt.figure(figsize=fig_dim) 
        if show_title plt.title("Time-domain representation of IR output pulse") end
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel(norm ? "I(z=end) (norm.)" : "I(z=end) (arb. units)")
        plt.plot(t*1e15, norm ? Maths.normbymax(It0_IR[:,end]) : It0_IR[:,end] , color="red", label=L"\tau_{FWHM}="*string(round(τ_UV*1e15, digits=1) )*"fs")
        plt.plot(t*1e15,norm ? Maths.normbymax(It0_IR_envelope[:,end]) : It0_IR_envelope[:,end], color="black", linestyle="--")

        if (read_IR & show_IR )==true  # overlay measured input pulse 
            plt.plot(t_in*1e15,norm ? Maths.normbymax(I_in)  : maximum(It0_envelope[:,1]).*Maths.normbymax(I_in), color="green")
        end

        plt.legend(loc="upper right")

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"time_domain_IR_output.pdf"))
            else     
                plt.savefig(joinpath(out_path,"time_domain_IR_output.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 11: plot spatiotemporal UV pulse at output
        rsym = Hankel.Rsymmetric(q)
        
        plt.figure(figsize=fig_dim)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.13)
        if show_title plt.title("Output UV pulse") end
        plt.xlabel("t (fs)")
        plt.ylabel("r (mm)")
        plt.ylim(minimum(rsym*1e3)/2, maximum(rsym*1e3)/2)
        plt.xlim(minimum(grid.t*1e15)/2, maximum(grid.t*1e15)/2)
        plt.pcolormesh(grid.t*1e15, rsym*1e3,norm ? Maths.normbymax(abs2.(Hankel.symmetric(Et_UV[:, :, end], q)')) : abs2.(Hankel.symmetric(Et_UV[:, :, end], q)'), cmap="inferno")
        #plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)"))
        #plt.yticks([])

        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"UV_pulse_output.pdf"))
		plt.savefig(joinpath(out_path,"UV_pulse_output.png"),dpi=1000)
		plt.savefig(joinpath(out_path,"UV_pulse_output.eps"))
            else     
                plt.savefig(joinpath(out_path,"UV_pulse_output.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 12: plot spatiotemporal IR pulse at output
        rsym = Hankel.Rsymmetric(q)
        
        plt.figure(figsize=fig_dim)
        plt.subplots_adjust(left=0.15, right=0.95, bottom=0.13)
        if show_title plt.title("Output IR pulse") end
        plt.xlabel("t (fs)")
        plt.ylabel("r (mm)")
        plt.ylim(minimum(rsym*1e3)/2, maximum(rsym*1e3)/2)
        plt.xlim(minimum(grid.t*1e15)/2, maximum(grid.t*1e15)/2)
        plt.pcolormesh(grid.t*1e15, rsym*1e3,norm ? Maths.normbymax(abs2.(Hankel.symmetric(Et_IR[:, :, end], q)')) : abs2.(Hankel.symmetric(Et_IR[:, :, end], q)'), cmap="inferno")
        plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)"))
        #plt.yticks([])

        plt.tight_layout()

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"IR_pulse_output.pdf"))
		plt.savefig(joinpath(out_path,"IR_pulse_output.png"),dpi=1000)
		plt.savefig(joinpath(out_path,"IR_pulse_output.eps"))
            else     
                plt.savefig(joinpath(out_path,"IR_pulse_output.png"),dpi=1000)
            end 
        end
        
        #+++++ PLOT 13: plot spatiotemporal UV pulse evolution 
        if read_ρ==true 
            z_vals_local =[L_total/2,L_total/2+L*0.5,L_total/2+L,L_total] # re-define zvals 
        else 
            z_vals_local = z_vals 
        end       
        idcs = [argmin(abs.(zout .- k)) for k in z_vals_local]   
        hide_cbar = true # if true: hide colour bars
        
        plt.figure(figsize=fig_dim) 
        if show_title plt.suptitle("Spatiotemporal evolution of UV pulse") end
        plt.subplots_adjust(hspace=0.6, wspace=0.55, bottom=0.13)

        for i in 1:4
            plt.subplot(2,2,i)
            plt.title("z="*string(round(z_vals_local[i]*1e3, digits=3))*"mm")
            plt.xlabel("t (fs)")
            plt.ylabel("r (mm)")
            plt.ylim(minimum(rsym*1e3)/2, maximum(rsym*1e3)/2)
            plt.xlim(minimum(grid.t*1e15)/2, maximum(grid.t*1e15)/2)
            plt.pcolormesh(grid.t*1e15, rsym*1e3,norm ? Maths.normbymax(abs2.(Hankel.symmetric(Et_UV[:, :, idcs[i]], q)')) : abs2.(Hankel.symmetric(Et_UV[:, :, idcs[i]], q)'))
            if hide_cbar==false  plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)")) end 
        end

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"UV_pulse_evolution.pdf"))
            else     
                plt.savefig(joinpath(out_path,"UV_pulse_evolution.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 14: plot spatiotemporal IR pulse evolution 
        plt.figure(figsize=fig_dim) 
        if show_title plt.suptitle("Spatiotemporal evolution of IR pulse") end
        plt.subplots_adjust(hspace=0.6, wspace=0.55, bottom=0.13)

        for i in 1:4
            plt.subplot(2,2,i)
            plt.title("z="*string(round(z_vals_local[i]*1e3, digits=3))*"mm")
            plt.xlabel("t (fs)")
            plt.ylabel("r (mm)")
            plt.ylim(minimum(rsym*1e3)/2, maximum(rsym*1e3)/2)
            plt.xlim(minimum(grid.t*1e15)/2, maximum(grid.t*1e15)/2)
            plt.pcolormesh(grid.t*1e15, rsym*1e3,norm ? Maths.normbymax(abs2.(Hankel.symmetric(Et_IR[:, :, idcs[i]], q)')) : abs2.(Hankel.symmetric(Et_IR[:, :, idcs[i]], q)'))
            if hide_cbar==false  plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)")) end
        end

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"IR_pulse_evolution.pdf"))
            else     
                plt.savefig(joinpath(out_path,"IR_pulse_evolution.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 15: UV spectral evolution
        c = [plt.get_cmap("viridis")(i) for i in range(0,1,length(z_vals_local))]
        
        if norm 
            max = 1
            for i=1:length(z_vals_local)
                nn = maximum(Iω0[λlowidx:λhighidx,idcs[i]])
                if (nn > max) max=nn end 
            end 
        end     
        
        maxmax = maximum(maximum(Iω0[λlowidx:λhighidx,:]))

        plt.figure(figsize=fig_dim)
        if show_title plt.title("UV beam spectral evolution") end
        for i in 1:length(z_vals_local)
            plt.plot(λ[λlowidx:λhighidx]*1e9,norm ? Iω0[λlowidx:λhighidx,idcs[i]]/max : Iω0[λlowidx:λhighidx,idcs[i]], label="z="*string(round(z_vals_local[i]*1e3, digits=3))*"mm", color=c[i])
        end
        plt.xlim(170, 340)
        plt.xlabel(L"\lambda"*" (nm)")
        plt.ylabel(norm ? "I (norm.)" : "I (arb. units)")

        plt.legend(loc="upper right")

        if save==true
            if use_pdf == true 
                plt.savefig(joinpath(out_path,"UV_spectral_evolution.pdf"))
            else     
                plt.savefig(joinpath(out_path,"UV_spectral_evolution.png"),dpi=1000)
            end 
        end

        #+++++ PLOT 16: spectral phase IR
        plt.figure(figsize=fig_dim)
        if show_title plt.title("IR spectral phase") end
         
        plt.plot(ω[ωlowIRidx:ωhighIRidx]*1e-15,ϕω0[ωlowIRidx:ωhighIRidx,1], color="grey", label="z="*string(round(zout[1]*1e3, digits=3))*"mm")
        plt.plot(ω[ωlowIRidx:ωhighIRidx]*1e-15,ϕω0[ωlowIRidx:ωhighIRidx,end], color="red", label="z="*string(round(zout[end]*1e3, digits=3))*"mm")
         
         plt.ylabel(L"\varphi"*" (rad)")
         plt.xlabel(L"\omega"*" "*L"(10^{15}"*"rad/s)")
         plt.legend(loc="upper right")

         if save==true
             if use_pdf == true 
                 plt.savefig(joinpath(out_path,"spectral_phase_IR.pdf"))
             else     
                 plt.savefig(joinpath(out_path,"spectral_phase_IR.png"),dpi=1000)
             end 
         end

        #+++++ PLOT 17: spectral phase UV output
        plt.figure(figsize=fig_dim)
        if show_title plt.title("UV spectral phase at output") end
         
	#plt.plot(λ[2:end]*1e9,norm ? Iω0[2:end,end]/maximum(Iω0[2:end,1]) : Iω0[2:end,end], label="z=$(round(zout[end]*1e3,digits=2))mm", color="red")
	#plt.plot(ω[ωlowUVidx:ωhighUVidx]*1e-15,Iω0[2:end,end], label="spectrum, z=$(round(zout[end]*1e3,digits=2))mm", color="gray")
	plt.plot(ω[ωlowUVidx:ωhighUVidx]*1e-15,ϕω0[ωlowUVidx:ωhighUVidx,end], label="spectral phase", color="orange")
         
        plt.ylabel(L"\varphi"*" (rad)")
        plt.xlabel(L"\omega"*" "*L"(10^{15}"*"rad/s)")

        plt.tight_layout()
 
         if save==true
             if use_pdf == true 
                 plt.savefig(joinpath(out_path,"spectral_phase_UV.pdf"))
             else     
                 plt.savefig(joinpath(out_path,"spectral_phase_UV.png"),dpi=1000)
             end 
         end 

        # +++++ PLOT 18:  Frequency-resolved UV output temporal profile 
        plt.figure(figsize=fig_dim)
        if show_title plt.title("Frequency-resolved UV output temporal profile") end
        
        function meshgrid(x, y)
            X = [i for i in x, j in 1:length(y)]
            Y = [j for i in 1:length(x), j in y]
            return X, Y
        end
        
        X, Y = meshgrid(λ[λlowidx:n:λhighidx]*1e9, t*1e15)    # n: bin size

        plt.contourf(X,Y, norm ? Maths.normbymax(I_ωt_UV) : I_ωt_UV, 10)
        plt.ylim(minimum(t)/2 *1e15, maximum(t)/2 *1e15)
        plt.colorbar(label=(norm ? "I (norm.)" : "I (arb. units)" ))
        plt.ylabel("t (fs)")
        plt.xlabel(L"\lambda"*"(nm)")

        plt.tight_layout()
        
        if save==true
            if use_pdf == true
                plt.savefig(joinpath(out_path,"UV_output_2d_map_time_frequency.pdf"))
            else 
                plt.savefig(joinpath(out_path,"UV_output_2d_map_time_frequency.png"),dpi=1000)
            end
        end  

    end    

    # ----------------- WRITE PARAMS & UV SPECTRUM TO FILE ------------------
    p_const ? dens_mod="const" : read_ρ ? dens_mod="coms" : dens_mod="grad"

    if save == true 
        open(joinpath(out_path,"params.txt"), "w") do file
            write(file, "gas     = "*string(gas)*"\n")
            if p_scan == false 
                write(file, "pres    = "*string(pres)*"\n")
            end    
            write(file, "p_ed    = "*string(p_ed)*"\n")
            write(file, "dens_mod = "*string(dens_mod)*"\n")
            write(file, "τ       = "*string(τ_input)*"\n")
            write(file, "λ0      = "*string(λ0)*"\n")
            write(file, "w0      = "*string(w0)*"\n")
            write(file, "CEP     = "*string(CEP)*"\n")
            write(file, "IRenergy  = "*string(IRenergy)*"\n")
            write(file, "L       = "*string(L)*"\n")
            write(file, "kerr    = "*string(kerr)*"\n")
            write(file, "ion     = "*string(ion)*"\n")
            write(file, "only_r0     = "*string(only_r0)*"\n")
            write(file, "propz   = "*string(propz)*"\n")
            write(file, "GVD     = "*string(ϕs[3])*"\n")
            write(file, "thickness= "*string(thickness)*"\n")
            write(file, "ion_mod = "*string(ion_model)*"\n")
            write(file, "N = "*string(N)*"\n")
            write(file, "R = "*string(R)*"\n")

            if p_scan == false
                write(file, "\n")
                write(file, "\n")
                write(file, "E_out   = "*string(UV_pulse_en[end])*"\n")
                write(file, "η       = "*string(η_THG)*"\n")
                write(file, "τ_UV    = "*string(τ_UV)*"\n")
                write(file, "z_peak  = "*string(z_peak)*"\n")
            end    
            
            if read_IR == true
                write(file, "\n")
                write(file, "# IR input beam read from: "*string(path_IR)*"\n")
            end

            if read_ρ == true
                write(file, "\n")
                write(file, "# Gas density distribution read from: "*string(path_ρ)*"\n")
            end

            if read_UV == true
                write(file, "\n")
                write(file, "# UV output beam read from: "*string(path_UV)*"\n")
            end
        end

        if p_scan == true 
            open(joinpath(out_path,string(pres)*" bar.dat"), "w") do file
                writedlm(file, zip(λ[λlowidx:λhighidx], Iω0[λlowidx:λhighidx,end]))
            end
        else     
            open(joinpath(out_path,"UV_spectrum.txt"), "w") do file
                writedlm(file, zip(λ[λlowidx:λhighidx], Iω0[λlowidx:λhighidx,end]))
            end
        
            if save_UV_temp==true 
                open(joinpath(out_path,"UV_temporal.txt"), "w") do file
                    writedlm(file, zip(t,It0_UV_envelope[:,end]))
                end    
            end

            if save_UV_phase==true 
                open(joinpath(out_path,"UV_phase.txt"), "w") do file
                    writedlm(file, zip(ω[ωlowUVidx:ωhighUVidx],ϕω0[ωlowUVidx:ωhighUVidx,end]))
                end    
            end
            if save_phase==true 
                open(joinpath(out_path,"phase_input.txt"), "w") do file
                    writedlm(file, zip(ω,ϕω0[:,1]))
                end
                open(joinpath(out_path,"spectrum_input.txt"), "w") do file
                    writedlm(file, zip(λ, Iω0[:,1]))
                end
                open(joinpath(out_path,"phase_output.txt"), "w") do file
                    writedlm(file, zip(ω,ϕω0[:,end]))
                end 
                open(joinpath(out_path,"spectrum_output.txt"), "w") do file
                    writedlm(file, zip(λ, Iω0[:,end]))
                end
            end
        end     
    end

    # ----------------- SHOW PLOTS ----------------------------
    if (show & !txt_only)==true 
        plt.show()
    end  

    # ----------------- RETURN OUTPUT ----------------------------
    if p_scan == true 
        return UV_pulse_en[end], η_THG, τ_UV, z_peak
    end
end    

# ----------------- RUN CODE ----------------------------

if p_scan == true 

    printstyled("Starting pressure scan on "*Dates.format(Dates.now(), "dd/mm/yyyy")*" at "*Dates.format(Dates.now(), "HH:MM:SS")*"\n", bold=true, color=:red, underline=true)

    out_path = joinpath("output", scan_dir)

    if isdir(out_path) == false 
        mkdir(out_path)
    end
        
    for i in range(1,length(pres_arr)) 
        printstyled("Starting run "*string(i)*"/"*string(length(pres_arr))*" at "*Dates.format(Dates.now(), "HH:MM:SS")*"\n", bold=true, color=:green)

        E_UV, η, τ_UV, z_peak = THG_main(pres_arr[i])

        open(joinpath(out_path,"energy_efficiency_time_zpeak.txt"), "a") do file
            writedlm(file, zip(pres_arr[i], E_UV, η, τ_UV, z_peak))
        end
    end 
    
    printstyled("Finished pressure scan on "*Dates.format(Dates.now(), "dd/mm/yyyy")*" at "*Dates.format(Dates.now(), "HH:MM:SS")*"\n", bold=true, color=:red, underline=true)

else 
    THG_main()
end     
