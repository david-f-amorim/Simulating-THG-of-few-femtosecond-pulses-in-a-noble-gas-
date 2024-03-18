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
txt_only = true             # if true, no plots are produced (only simulation parameters, UV spectrum and phases are written to file)

read_IR = true               # if true: read input IR pulse [time domain] from file; if false: use Gaussian approximation 
read_ρ  = false              # if true: read gas density profile from file; if false: use pressure gradient approximation 
read_UV = false              # if true: read in measured UV spectrum and overlay on simulated results   
only_r0 = false		     # if false: E-field is integrated over r for time domain analysis, intensity analysis, ...

IR_spec = false              # if true: read measured input IR spectrum from file and overlay 
show_IR = false              # if true and "read_IR" is true: overlay measured time-domain input pulse on plots (barely ever relevant; used to check that transform to frequency domain works)
IR_spec_exp = false          # if true: read input IR spectrometer spectrum from file and overlay (barely ever relevant; used to check that FROG spectrum is accurate)

save_UV_temp = true          # if true: write time-domain output UV pulse to file (has no effect for pressure scans)
save_UV_phase= true
save_phase = true	    # if true: write full spectral phase at input and output to two files

# ------------------ SET PHYSICAL PARAMETERS ------------------------

gas = :Ne           # gas type in cell 
pres = 5.0          # central gas pressure [bar] (for single run; not relevant for pressure scans)
p_ed = 1e-3         # edge gas pressure [bar] 
p_const = true     # if true: set constant pressure profile P==(pres,pres,pres) ; if false: set simple gradient: P==(p_ed, pres, p_ed); (only relevant when read_ρ==false)
τ = 5e-15           # FWHM pulse duration of IR pulse [s] (only relevant when read_IR==false)
λ0 = 800e-9         # central wavelength of IR pulse [m]
w0 = 65e-6          # beam waist of IR pulse [m]
CEP = 0.0           # carrier-envelope phase of IR pulse (ϕ0) [rad] (only relevant when read_IR==true)                                 
IRenergy = 400e-6    # IR pulse energy [J] (related to beam power via 1kHz repetition rate)                                                           
L = 5e-3            # propagation distance (cell length) [m]

propz = -L/2        # propagation distance from the waist [m], i.e. beam focus position  (NOTE: always specified in a coordinate system where the cell starts at 0 and ends at L!)
z_vals =L .* [1.0/4, 1.0/3, 3.0/4, 1.0]     # points along the cell at which to investigate beam evolution [m] (NOTE: always specified in a coordinate system where the cell starts at 0 and ends at L!)

λ_lims = (100e-9, 1000e-9)       # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (100e-9, 360e-9)     # wavelength limits of UV region of interest [m,m]
λ_rangeIR = (600e-9, 1000e-9)    # wavelength limits of IR region of interest [m,m]

material = :SiO2     # material of the optics the IR beam propagates through before THG (for chirp compensation)
thickness= 0         # thickness of the material [m] (for chirp compensation)

ϕs = [0,0,0,0]     # Taylor-series coefficients of initial spectral phase  [s^n] (used to introduce additional chirp)

ion = false         # if true: enable ionisation response, if false: disable ionisation 
ion_model="PPT"    # set to "ADK" or "PPT" (has no effect if ion==false); "ADK" is less accurate at low intensities but faster; "PPT" may crash at very high intensities

pres_arr = range(start= 0.25, stop= 2.0, step= 0.25)  # pressure range (only relevant for pressure scans)  [bar]


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

        file_ρ    = "dens_$(pres)bar.dat"                     # name of density profile data file 
        path_ρ    = joinpath(in_dir, file_ρ)                  # sys. path to density profile data file 
    
        z_in = readdlm(path_ρ,' ', Float64, '\n')[:,1]        # read in spatial points 
        ρ_in = readdlm(path_ρ,' ', Float64, '\n')[:,2]        # read in density data 
        
        dens = Maths.CSpline(z_in, ρ_in)                      # interpolate density function 
        γ = PhysData.sellmeier_gas(gas)                       # get polarisability from Sellmeier expansion
        
        coren(ω; z) = sqrt(1 + γ(wlfreq(ω)*1e6)*dens(z))      # calculate refractive index along the cell   

        L_total = maximum(z_in)                               # get total propagation distance 
        
    else   
        Z= (0, L/2, L)                                        # define points for pressure gradient
        p_const ? P=(pres,pres,pres) : P= (p_ed, pres, p_ed)  # define values for pressure gradient (see definition "p_const" above)
        (coren,dens)=Capillary.gradient(gas,Z,P)              # gives n(ω; z) and ρ(z) from pressure profile   
    end     

    # ----------------- SET SIMULATION GRID ----------------------------

    R = 200e-6            # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  
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
    ϕω0 = angle.(Er0)                                # spectral phase of the complete field across all r 
    ϕωr = angle.(Erout)                              # spectral phase of the complete field at r≠0

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

        

    end    

    # ----------------- WRITE PARAMS & UV SPECTRUM TO FILE ------------------
    read_ρ ? dens_mod="coms" : dens_mod="grad"

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
        end     # open and write "params.txt" file

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
                open(joinpath(out_path,"spectrum.txt"), "w") do file
                    writedlm(file, zip(λ[λlowidx:λhighidx], Iω0[λlowidx:λhighidx,end]))
                end
                open(joinpath(out_path,"phase_output.txt"), "w") do file
                    writedlm(file, zip(ω,ϕω0[:,end]))
                end    
            end
        end    # if/else p_scan==true 
    end     # if save==true

    # ----------------- SHOW PLOTS ----------------------------
    if (show & !txt_only)==true 
        plt.show()
    end  

    # ----------------- RETURN OUTPUT ----------------------------
    if p_scan == true 
        return UV_pulse_en[end], η_THG, τ_UV, z_peak
    end
end  # THG_main() function  

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
