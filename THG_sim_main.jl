using  Luna
import Luna.PhysData: wlfreq    
import FFTW                    
import Luna: Hankel            
import Dates                   
using  DelimitedFiles         

# ----------------- QUICK SETTINGS -------------------------------
p_scan = false                # if true, run is counted as part of a pressure scan 

save = true                  # if true, saves output plots, run parameters and UV spectrum
show = true                  # if true, opens plots in GUI after run 
txt_only = false             # if true, no plots are produced 

read_IR = true               # if true: read input IR pulse from file; if false: use Gaussian approximation 
read_ρ  = false              # if true: read gas density profile from file; if false: use pressure gradient approximation 
read_UV = false              # if true: overlay measured UV output on simulated results   

IR_spec = false              # if true: read input IR spectrum from file and overlay 
show_IR = false              # if true and "read_IR" is true: overlay measured input pulse on plots

# ------------------ SET MEASURED PARAMETERS ------------------------

gas = :Ar          # gas
pres = 0.4          # central gas pressure [bar]   (ignored if p_scan is true)
p_ed = 1e-3         # edge gas pressure [bar]
p_const = false     # if true: set constant pressure profile P==(pres,pres,pres) ; if false: set simple gradient: P==(p_ed, pres, p_ed)
τ = 5e-15           # FWHM pulse duration [s] (only relevant when temporal beam profile is approximated as Gaussian)
λ0 = 800e-9         # central wavelength [m]
w0 = 65e-6          # beam waist [m]
ϕ = 0.0             # carrier-envelope offset (CEO) phase [rad]                                    -> can this be extracted from data?
energy = 150e-6     # pulse energy [J]                                                             -> multiply by 1kHz (?) repetition rate for beam power
L = 3e-3            # propagation distance (cell length) [m]

λ_lims = (100e-9, 1000e-9)       # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (100e-9, 360e-9)     # wavelength limits of UV region of interest [m,m]
λ_rangeIR = (600e-9, 1000e-9)    # wavelength limits of IR region of interest [m,m]

# ------------------ SET NONLINEAR PARAMS ------------------------

kerr = "f"         # set nonlinear Kerr effect: 
                    #       must be "ff" [Kerr_field + Kerr_field_nothg], "ef" [Kerr_env + Kerr_field_nothg], 
                    #       "f" [Kerr_field], "fe" [Kerr_field + Kerr_env] , or "e" [Kerr_env]
ion = true          # if true: enable ionisation response

# ---------------- SET SCAN PARAMETERS ------------------------

scan_dir = "scan_"*string(energy*1e6)*"mW_"*string(gas)*"_"*string(round(ϕ; digits=3))*"rad_"*string(kerr)*"_"*string(ion ? "ion" : "no-ion")  # name of scan output directory 

pres_arr = range(start= 0.1, stop= 5.1, step= 0.1)  # pressure range [bar]

# ----------------- INPUT HANDLING -------------------------------

# NOTE: see "./file_prepare.py" for relevant file content specifications

in_dir    = "input"                          # directory of input files

file_IR   = "IRpulse.dat"                    # name of IR input pulse file 
path_IR   = joinpath(in_dir, file_IR)        # sys. path to IR input pulse file 

file_ρ    = "dens.dat"                       # name of density profile data file 
path_ρ    = joinpath(in_dir, file_ρ)         # sys. path to density profile data file 

file_UV   = "UVpulse.dat"                    # name of UV output pulse file 
path_UV   = joinpath(in_dir, file_UV)        # sys. path to UV output pulse file 

file_IR_spec = "IRspec.dat"                   # name of IR input spectrum file 
path_IR_spec = joinpath(in_dir, file_IR_spec) # sys. path to IR input spectrum file 

# ----------------- MAKE ARRANGEMENTS FOR PRESSURE SCAN -----------

if p_scan == true  
    txt_only = true 
    save     = true 
end    

# ----------------- DEFINE MAIN FUNCTION -----------

function THG_main(pres=pres) 

    # ----------------- SET PRESSURE PROFILE ---------------------------

    if read_ρ == true 
        z_in = readdlm(path_ρ,' ', Float64, '\n')[:,1]        # read in spatial points 
        ρ_in = readdlm(path_ρ,' ', Float64, '\n')[:,2]        # read in density data 
        
        dens = Maths.CSpline(z_in, ρ_in)                      # interpolate density function 
        γ = PhysData.sellmeier_gas(gas)                       # get polarisability from Sellmeier expansion
        
        coren(ω; z) = sqrt(1 + γ(wlfreq(ω)*1e6)*dens(z))      # calculate refractive index along the cell   

    else   
        Z= (0, L/2, L)                                        # define points for pressure gradient
        p_const ? P=(pres,pres,pres) : P= (p_ed, pres, p_ed)  # define values for pressure gradient (see "p_const" above)
        (coren,dens)=Capillary.gradient(gas,Z,P)              # gives n(ω; z) and ρ(z) from pressure profile   
    end     

    

    # ----------------- SET SIMULATION GRID ----------------------------

    R = 250e-6            # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  
    N = 1024              # sample size for Hankel tansform grid     [1]                         ->> WHY THIS VALUE?  WOULD CHOOSING N≫ IMPROVE ACCURACY?
    trange = 0.05e-12     # total extent of time window required [s]                             


    grid = Grid.RealGrid(L, λ0, λ_lims ,trange)   # set up time grid
    q = Hankel.QDHT(R, N, dim=2)                  # set up discrete Hankel transform matrix                        

    energyfun, _ = Fields.energyfuncs(grid, q)    # "energyfun" gives total energy in a field E(t)    

    # ----------------- SET NONLINEAR EFFECTS ----------------------------

    ionpot = PhysData.ionisation_potential(gas)                 # set gas ionisation potential   
    ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)        # set gas ionisation rate 

    linop = LinearOps.make_linop(grid, q, coren)                # generate linear operator for pulse-propagation equation  
    normfun = NonlinearRHS.norm_radial(grid, q, coren)          # generate normalisation function for radial symmetry 

    # * * * SET KERR EFFECT AND/OR PLASMA FORMATION 
    if (kerr=="ff") & (ion == true)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas),length(grid.to))
                )
    elseif (kerr=="ef") & (ion == true)
        responses = (Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas),length(grid.to))
                )
    elseif (kerr=="f") & (ion == true)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                )
    elseif (kerr=="e") & (ion == true)
        responses = (Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                )
    elseif (kerr=="fe") & (ion == true)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot),
                Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                )   
    elseif (kerr=="ff") & (ion == false)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas),length(grid.to))
                )
    elseif (kerr=="ef") & (ion == false)
        responses = (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                Nonlinear.Kerr_field_nothg(PhysData.γ3_gas(gas),length(grid.to))
                )    
    elseif (kerr=="f") & (ion == false)
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                )    
    elseif (kerr=="e") & (ion == false)
        responses = (Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                )    
    elseif (kerr=="fe") & (ion == false) 
        responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),     
                Nonlinear.Kerr_env(PhysData.γ3_gas(gas)),
                )
    end    

    # ----------------- SET INPUT FIELD ----------------------------                                       

    if read_IR == true 

        t_in = readdlm(path_IR,' ', Float64, '\n')[:,1]                   # read in timing data
        I_in = readdlm(path_IR,' ', Float64, '\n')[:,2]                   # read in intensity data (time domain)

        beam_spline    = Maths.CSpline(t_in, I_in)                           # interpolates temporal beam intensity envelope  
        beam_spline_fs = Maths.CSpline(t_in*1e15, I_in)                      # interpolates temporal beam intensity envelope with t in fs

        function beam_profile(t,r::AbstractVector)
        beam_spline.(t) .* Maths.gauss.(r, w0/2)'                      # models spatial beam profile as Gaussian  
        end

        inputs = Fields.SpatioTemporalField(λ0, energy, ϕ, 0.0,           # models input beam based off measured data
        (t, r) -> beam_profile(t, r), -L/2)  

    else 
        inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ,                   # models input beam as Gaussian 
                energy=energy, w0=w0, propz=-L/2)
    end            

    # ----------------- RUN SIMULATION ----------------------------

    Eω, transform, FT = Luna.setup(grid, q, dens, normfun, responses, inputs)  # set-up propagation 
    output = Output.MemoryOutput(0, grid.zmax, 201)                            # configure output      ->> WHY THESE VALUES ???!!!!
    Luna.run(Eω, grid, linop, transform, FT, output)                           # run simulation  

    # ----------------- PROCESS RESULTS ----------------------------

    # * * * EXTRACT GRID PARAMS:  
    ω = grid.ω                   # sampled angular frequency values [rad/s]
    f = ω / 2π                   # sampled linear frequency values [Hz]
    λ = ω .\ (2π * PhysData.c)   # sampled wavelengths [m]                 ->> CAREFUL: λ[1]=Inf !
    t = grid.t                   # sampled points in time [s]                                          
    zout = output.data["z"]      # sampled points along z [m] 

    # * * * DEFINE USEFUL GRID INDICES:
    ω0idx         = argmin(abs.(ω .- 2π*PhysData.c/λ0))                   # index X such that ω[X]=ω0 (fundamental)
    ωTHidx        = argmin(abs.(ω .- 2π*PhysData.c/(λ0/3)))               # index X such that ω[X]=3ω0 (TH)
    ωhighUVidx    = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeUV[1]))         # index X such that ω[X]=ω_maxUV 
    ωlowUVidx     = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeUV[2]))         # index X such that ω[X]=ω_minUV
    ωhighIRidx    = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeIR[1]))         # index X such that ω[X]=ω_maxIR 
    ωlowIRidx     = argmin(abs.(ω .- 2π*PhysData.c/λ_rangeIR[2]))         # index X such that ω[X]=ω_minIR


    λhighidx    = argmin(abs.(λ .- λ_rangeUV[1]))         # index X such that λ[X]=λ_maxUV 
    λlowidx     = argmin(abs.(λ .-λ_rangeUV[2]))          # index X such that λ[X]=λ_minUV

    # * * * EXTRACT FIELD AMPLITUDES:
    Eout = output.data["Eω"]                         # Hankel-transformed amplitude in frequency domain  
    Erout = (q \ Eout)                               # Real-space amplitude in frequency domain at r ≠ 0   ->> "\" represents inverse Hankel transform 
    Er0 = dropdims(Hankel.onaxis(Eout, q), dims=2)   # Real-space amplitude in frequency domain at r=0 

    Etout = FFTW.irfft(Erout, length(t), 1)     # time-domain real field amplitude at r≠0
    Et0 = FFTW.irfft(Er0, length(t),1)          # time-domain real field amplitude at r=0

    # * * * CONVERT TO INTENSITIES:
    Iωr = abs2.(Erout)                               #  intensity at r≠0 in frequency domain  [arbitrary units]
    Iω0 = abs2.(Er0)                                 #  intensity at r=0 in frequency domain  [arbitrary units]
    It0 = abs2.(Et0)                                 #  intensity at r=0 in time domain [arbitrary units]

    # * * * FILTER FOR UV FIELD (r≠0):
    filter=FFTW.rfft(Etout, 1)        # set up filter array 
    filter[1:ωlowUVidx,:,:].=0;       # filters out ω < ω_min
    filter[ωhighUVidx:end,:,:].=0;    # filters out ω > ω_max 

    Etout_UV=FFTW.irfft(filter, length(t), 1)    # time-domain real field amplitude of UV pulse 

    if txt_only == false 
        # * * * FILTER FOR UV FIELD (r=0)
        filter_onaxis = FFTW.rfft(Et0, 1)          # set up filter array
        filter_onaxis[1:ωlowUVidx,:].=0;           # filters out ω < ω_min
        filter_onaxis[ωhighUVidx:end,:].=0;        # filters out ω > ω_max  

        Et0_UV =FFTW.irfft(filter_onaxis, length(t), 1)    # time-domain real field amplitude of UV pulse at r=0
        It0_UV = abs2.(Et0_UV)                           # intensity of on-axis UV pulse at r=0

        # * * * FILTER FOR IR FIELD (r=0)
        filter_onaxis_IR = FFTW.rfft(Et0, 1)          # set up filter array
        filter_onaxis_IR[1:ωlowIRidx,:].=0;           # filters out ω < ω_min
        filter_onaxis_IR[ωhighIRidx:end,:].=0;        # filters out ω > ω_max  

        Et0_IR =FFTW.irfft(filter_onaxis_IR, length(t), 1)    # time-domain real field amplitude of IR pulse at r=0
    
        # * * * EXTRACT INTENSITY ENVELOPES 
        It0_envelope = abs2.(Maths.hilbert(Et0))          # envelope modulating It0
        It0_UV_envelope = abs2.(Maths.hilbert(Et0_UV))    # envelope modulating It0_UV
        It0_IR_envelope = abs2.(Maths.hilbert(Et0_IR))    # envelope modulating It0_IR 

        # * * * CALCULATE PULSE DURATIONS
        τ_input = Maths.fwhm(t, It0_envelope[:,1])             # FWHM input pulse duration     [s]
        τ_UV    = Maths.fwhm(t, It0_UV_envelope[:,end])        # FWHM UV output pulse duration [s]
    end

    # * * * CALCULATE PULSE ENERGIES:
    tot_pulse_en = zeros(length(zout))                   # set up array for total pulse energy 
    UV_pulse_en  = zeros(length(zout))                   # set up array for UV pulse energy 

    for i = 1:size(Etout, 3)
        tot_pulse_en[i] = energyfun(Etout[:, :, i])         # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
        UV_pulse_en[i]  = energyfun(Etout_UV[:, :, i])      # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
    end     

    η_THG = UV_pulse_en[end]/tot_pulse_en[1]            # THG efficiency: initial total pulse energy / final UV pulse energy 

    # * * * PROCESS MEASURED DATA FROM FILES 
    if (IR_spec == true) & (txt_only == false)
        λ_IR_spec = readdlm(path_IR_spec,' ', Float64, '\n')[:,1]             # read in IR input wavelength data 
        I_IR_spec = readdlm(path_IR_spec,' ', Float64, '\n')[:,2]             # read in IR input spectral data 
    end

    if (read_UV == true) & (txt_only == false) 
        λ_in    = readdlm(path_UV,' ', Float64, '\n')[:,1]                      # read in UV wavelength data
        I_in_UV = readdlm(path_UV,' ', Float64, '\n')[:,2]                      # read in UV intensity data 
    end    

    # ----------------- OUTPUT HANDLING -------------------------

    if isdir("output")== false   # create general directory to store output
        mkdir("output")
    end    

    if (save & !p_scan) == true              # create specific directory to store run output 
        out_path = joinpath("output","run_"*Dates.format(Dates.now(), "yyyy_mm_dd__HH_MM"))
        mkdir(out_path)
    elseif p_scan == true  
        out_path = joinpath("output", scan_dir)
    end 

    # ----------------- PLOT RESULTS ----------------------------
    if txt_only == false 

        @eval import PyPlot:pygui, plt
        close("all")
        pygui(true)

        #+++++ PLOT 1:  fundamental and third harmonic intensities as functions of z and r≠0
        plt.figure(figsize=[7.04, 5.28])
        plt.suptitle("Off-axis intensity of fundamental and third harmonic")
        plt.subplots_adjust(left=0.125, bottom=0.11, right=0.992, top=0.88, hspace=0.5)

        plt.subplot(2,1,1)
        plt.pcolormesh(zout*1e3, q.r*1e3, Iωr[ω0idx, :, :])
        plt.colorbar(label="arb. units")
        plt.ylabel("r (mm)")
        plt.xlabel("z (mm)")
        plt.title("I(r,z; λ=$(round(Int,λ0*1e9))nm)")    

        plt.subplot(2,1,2)
        plt.pcolormesh(zout*1e3, q.r*1e3,Iωr[ωTHidx, :, :])
        plt.colorbar(label="arb. units")
        plt.xlabel("z (mm)")
        plt.ylabel("r (mm)")
        plt.title("I(r,z; λ=$(round(Int,λ0/3*1e9))nm)")

        if save==true
            plt.savefig(joinpath(out_path,"off-axis_intensity.png"),dpi=1000)
        end    
            
        #+++++ PLOT 2:  fundamental and third harmonic intensities as functions of z at r=0
        plt.figure(figsize=[7.04, 5.28])
        plt.suptitle("On-axis intensity of fundamental and third harmonic")
        plt.subplots_adjust(hspace=0.5)

        plt.subplot(2,1,1)
        plt.plot(zout*1e3,  Iω0[ω0idx,  :], color="red")
        plt.title("λ=$(round(Int,λ0*1e9))nm")
        plt.ylabel("I(r=0) (arb. units)")
        plt.xlabel("z (mm)")
        plt.ticklabel_format(axis="y", style="scientific", scilimits=(0,0))

        plt.subplot(2,1,2)
        plt.plot(zout*1e3,  Iω0[ωTHidx,  :], color="red")
        plt.title("λ=$(round(Int,λ0/3*1e9))nm")
        plt.xlabel("z (mm)")
        plt.ylabel("I(r=0) (arb. units)")
        plt.ticklabel_format(axis="y", style="scientific", scilimits=(0,0))

        if save==true
            plt.savefig(joinpath(out_path,"on-axis_intensity.png"),dpi=1000)
        end 

        #+++++ PLOT 3: gas number density and effective susceptibility along the cell 
        plt.figure(figsize=[7.04, 5.28])
        plt.suptitle("Initial gas properties")
        plt.subplots_adjust(hspace=0.4)

        plt.subplot(2,1,1)
        plt.plot(zout*1e3, [dens(i) for i in zout] ./ PhysData.N_A, label="central pressure: P₀=$(pres) bar", color="red")
        plt.ylabel("ρ (mol/m³)")
        plt.xlabel("z (mm)")
        plt.legend()

        plt.subplot(2,1,2)

        χ0  = [coren(ω[ω0idx],z=i)^2-1 for i in zout]
        χTH = [coren(ω[ωTHidx],z=i)^2-1 for i in zout]
        
        plt.xlabel("z (mm)")
        plt.ylabel("Effective linear χₑ")
        plt.plot(zout*1e3, χ0, label="λ=$(round(Int,λ0*1e9))nm", color="grey")
        plt.plot(zout*1e3, χTH, label="λ=$(round(Int,λ0/3*1e9))nm", color="red")
        plt.ticklabel_format(axis="y", style="scientific", scilimits=(0,0))
        plt.legend()

        if save==true
            plt.savefig(joinpath(out_path,"density_and_susceptibility.png"),dpi=1000)
        end 

        #+++++ PLOT 4:  linear on-axis spectrum I(λ) at z=0 and z=L 
        plt.figure(figsize=[7.04, 5.28])
        plt.title("Linear on-axis beam spectrum")
        plt.plot(λ[2:end]*1e9, Iω0[2:end,1], label="z=0mm", color="grey")
        plt.plot(λ[2:end]*1e9, Iω0[2:end,end], label="z=$(L*1e3)mm", color="red")
        plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
        plt.xlabel("λ (nm)")
        plt.ylabel("I(r=0, λ) (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9, maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9, maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec), color="green", ls="--") 
        end    

        plt.legend()

        if save==true
            plt.savefig(joinpath(out_path,"full_on-axis_spectrum.png"),dpi=1000)
        end 

        #+++++ PLOT 5:  UV only linear on-axis spectrum I(λ) at z=0 and z=L 
        plt.figure(figsize=[7.04, 5.28])
        plt.title("Linear on-axis UV spectrum")
        plt.plot(λ[λlowidx:λhighidx]*1e9, Iω0[λlowidx:λhighidx,1], label="z=0mm", color="grey")
        plt.plot(λ[λlowidx:λhighidx]*1e9, Iω0[λlowidx:λhighidx,end], label="z=$(L*1e3)mm", color="red")
        plt.xlim(λ_rangeUV[1]*1e9, λ_rangeUV[2]*1e9)
        plt.xlabel("λ (nm)")
        plt.ylabel("I(r=0, λ) (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9, maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9, maximum(Iω0[λlowidx:λhighidx,1]).*Maths.normbymax(I_IR_spec), color="green", ls="--") 
        end

        plt.legend()

        if save==true
            plt.savefig(joinpath(out_path,"UV_on-axis_spectrum.png"),dpi=1000)
        end 

        #+++++ PLOT 6:  log. on-axis spectrum I(λ) at z=0 and z=L 
        Iω0log = log10.(Iω0)

        plt.figure(figsize=[7.04, 5.28])
        plt.title("Logarithmic on-axis beam spectrum")
        plt.plot(λ[2:end]*1e9, Iω0log[2:end,1], label="z=0mm",  color="grey")
        plt.plot(λ[2:end]*1e9, Iω0log[2:end,end], label="z=$(L*1e3)mm", color="red")
        plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
        plt.xlabel("λ (nm)")
        plt.ylabel("I(r=0, λ) (arb. units)")

        if read_UV == true  # overlay measured UV output spectrum 
            plt.plot(λ_in*1e9, log10.(abs.(maximum(Iω0[λlowidx:λhighidx,end]).*Maths.normbymax(I_in_UV))), color="purple", label="UV data (rescaled)")
        end

        if IR_spec == true 
            plt.plot(λ_IR_spec*1e9, log10.(abs.(maximum(Iω0[2:end,1]).*Maths.normbymax(I_IR_spec))), color="green", ls="--") 
        end

        plt.legend()


        if save==true
            plt.savefig(joinpath(out_path,"on-axis_spectrum_log.png"),dpi=1000)
        end 

        #+++++ PLOT 7:  pulse energies and efficiency 
        plt.figure(figsize=[7.04, 5.28])
        plt.suptitle("THG conversion efficiency: η="*string(round(η_THG, digits=4) *100)*"%")
        plt.subplots_adjust(hspace=0.5)

        plt.subplot(2,1,1)
        plt.plot(zout.*1e3, tot_pulse_en.*1e6, label="ΔE=-$(round(Int64, tot_pulse_en[1]*1e6-tot_pulse_en[end]*1e6))μJ", color="red")
        plt.xlabel("z (mm)")
        plt.ylabel("E (μJ)")
        plt.title("Total pulse energy")
        plt.legend()

        plt.subplot(2,1,2)
        plt.plot(zout.*1e3, UV_pulse_en.*1e9, label="ΔE=+$(round(Int64, UV_pulse_en[end]*1e9))nJ", color="red")
        plt.xlabel("z (mm)")
        plt.ylabel("E (nJ)")
        plt.title("UV pulse energy")
        plt.legend()

        if save==true
            plt.savefig(joinpath(out_path,"pulse_energies.png"),dpi=1000)
        end

        #+++++ PLOT 8: on-axis frequency evolution 
        plt.figure(figsize=[7.04, 5.28])
        plt.suptitle("On-axis frequency evolution")
        plt.pcolormesh(zout*1e3, f*1e-15, log10.(Maths.normbymax(Iω0)))   
        plt.colorbar(label="arb. units, normed")
        plt.clim(0, -6)  
        plt.xlabel("z (mm)")
        plt.ylabel("f (PHz)")
        plt.title("log. I(r=0, ω)")

        if save==true
            plt.savefig(joinpath(out_path,"on-axis_frequency_evolution.png"),dpi=1000)
        end

        #+++++ PLOT 9: time-domain plot of input pulse 
        plt.figure(figsize=[7.04, 5.28]) 
        plt.title("Time-domain representation of on-axis input pulse")
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel("I(t; r=0, z=0) (arb. units)")
        plt.plot(t*1e15, It0[:,1] , color="red", label="FWHM pulse duration: τ="*string(round(τ_input*1e15, digits=1) )*"fs")
        plt.plot(t*1e15,It0_envelope[:,1], color="black", ls="--")

        if (read_IR & show_IR )==true  # overlay measured input pulse 
            plt.plot(t_in*1e15,maximum(It0_envelope[:,1]).*Maths.normbymax(I_in), color="green")
        end

        plt.legend(loc="upper right")

        if save==true
            plt.savefig(joinpath(out_path,"time_domain_input.png"),dpi=1000)
        end

        #+++++ PLOT 10: time-domain plot of UV output pulse 
        plt.figure(figsize=[7.04, 5.28]) 
        plt.title("Time-domain representation of on-axis UV output pulse")
        plt.xlabel("t (fs)")
        plt.ylabel("I(t; r=0, z=L) (arb. units)")
        plt.plot(t*1e15, It0_UV[:,end] , color="red", label="FWHM pulse duration: τ="*string(round(τ_UV*1e15, digits=1) )*"fs")
        plt.plot(t*1e15,It0_UV_envelope[:,end], color="black", linestyle="--")

        plt.legend(loc="upper right")

        if save==true
            plt.savefig(joinpath(out_path,"time_domain_UV_output.png"),dpi=1000)
        end

        #+++++ PLOT 11: time-domain pulse evolution 
        plt.figure(figsize=[7.04, 5.28]) 
        plt.title("Evolution of total time-domain pulse")
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel("I(t; r=0) (arb. units)")

        c = [plt.get_cmap("viridis")(i) for i in range(0,1,5)]
        
        plt.plot(t*1e15,It0_envelope[:,1], label="z=0.0mm", color=c[5])
        plt.plot(t*1e15,It0_envelope[:,Int(round(length(zout)/4))], label="z=$(1/4*L*1e3)mm", color=c[4])
        plt.plot(t*1e15,It0_envelope[:,Int(2*round(length(zout)/4))], label="z="*string(round(2/4*L*1e3, digits=3))*"mm", color=c[3])
        plt.plot(t*1e15,It0_envelope[:,Int(3*round(length(zout)/4))], label="z=$(3/4*L*1e3)mm", color=c[2])
        plt.plot(t*1e15,It0_envelope[:,end], label="z=$(L*1e3)mm", color=c[1])

        plt.legend(loc="upper right")

        if save==true
            plt.savefig(joinpath(out_path,"pulse_evolution.png"),dpi=1000)
        end

        #+++++ PLOT 12: time-domain UV pulse evolution 
        plt.figure(figsize=[7.04, 5.28]) 
        plt.title("Evolution of UV time-domain pulse")
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel("I(t; r=0) (arb. units)")
        
        plt.plot(t*1e15,It0_UV_envelope[:,1], label="z=0.0mm", color=c[5])
        plt.plot(t*1e15,It0_UV_envelope[:,Int(round(length(zout)/4))], label="z=$(1/4*L*1e3)mm", color=c[4])
        plt.plot(t*1e15,It0_UV_envelope[:,Int(2*round(length(zout)/4))], label="z="*string(round(2/4*L*1e3, digits=3))*"mm", color=c[3])
        plt.plot(t*1e15,It0_UV_envelope[:,Int(3*round(length(zout)/4))], label="z=$(3/4*L*1e3)mm", color=c[2])
        plt.plot(t*1e15,It0_UV_envelope[:,end], label="z=$(L*1e3)mm", color=c[1])

        plt.legend(loc="upper right")

        if save==true
            plt.savefig(joinpath(out_path,"UV_pulse_evolution.png"),dpi=1000)
        end

        #+++++ PLOT 13: time-domain UV pulse evolution 
        plt.figure(figsize=[7.04, 5.28]) 
        plt.title("Evolution of IR time-domain pulse")
        plt.xlabel("t (fs)")
        plt.xlim(minimum(t)*1e15, maximum(t)*1e15)
        plt.ylabel("I(t; r=0) (arb. units)")
        
        plt.plot(t*1e15,It0_IR_envelope[:,1], label="z=0.0mm", color=c[5])
        plt.plot(t*1e15,It0_IR_envelope[:,Int(round(length(zout)/4))], label="z=$(1/4*L*1e3)mm", color=c[4])
        plt.plot(t*1e15,It0_IR_envelope[:,Int(2*round(length(zout)/4))], label="z="*string(round(2/4*L*1e3, digits=3))*"mm", color=c[3])
        plt.plot(t*1e15,It0_IR_envelope[:,Int(3*round(length(zout)/4))], label="z=$(3/4*L*1e3)mm", color=c[2])
        plt.plot(t*1e15,It0_IR_envelope[:,end], label="z=$(L*1e3)mm", color=c[1])

        plt.legend(loc="upper right")

        if save==true
            plt.savefig(joinpath(out_path,"IR_pulse_evolution.png"),dpi=1000)
        end
    end    
    # ----------------- WRITE PARAMS & UV SPECTRUM TO FILE ------------------

    if save == true 
        open(joinpath(out_path,"params.txt"), "w") do file
            write(file, "gas     = "*string(gas)*"\n")
            if p_scan == false 
                write(file, "pres    = "*string(pres)*"\n")
            end    
            write(file, "p_ed    = "*string(p_ed)*"\n")
            write(file, "p_const = "*string(p_const)*"\n")
            write(file, "τ       = "*string(τ)*"\n")
            write(file, "λ0      = "*string(λ0)*"\n")
            write(file, "w0      = "*string(w0)*"\n")
            write(file, "ϕ       = "*string(ϕ)*"\n")
            write(file, "energy  = "*string(energy)*"\n")
            write(file, "L       = "*string(L)*"\n")
            write(file, "kerr    = "*string(kerr)*"\n")
            write(file, "ion     = "*string(ion)*"\n")
            
            if p_scan == false
                write(file, "\n")
                write(file, "\n")
                write(file, "E_out   = "*string(UV_pulse_en[end])*"\n")
                write(file, "η       = "*string(η_THG)*"\n")
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
        end     
    end

    # ----------------- SHOW PLOTS ----------------------------
    if (show & !txt_only)==true 
        plt.show()
    end  

    # ----------------- RETURN OUTPUT ----------------------------
    if p_scan == true 
        return UV_pulse_en[end], η_THG
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

        E_UV, η = THG_main(pres_arr[i])

        open(joinpath(out_path,"energy_efficiency.txt"), "a") do file
            writedlm(file, zip(pres_arr[i], E_UV, η))
        end
    end 
    
    printstyled("Finished pressure scan on "*Dates.format(Dates.now(), "dd/mm/yyyy")*" at "*Dates.format(Dates.now(), "HH:MM:SS")*"\n", bold=true, color=:red, underline=true)

else 
    THG_main()
end     