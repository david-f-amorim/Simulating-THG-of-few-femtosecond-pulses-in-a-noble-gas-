using  Luna
import Luna.PhysData: wlfreq
import FFTW
import Luna: Hankel
import NumericalIntegration: integrate, SimpsonEven
import Dates
using  DelimitedFiles

# ----------------- QUICK SETTINGS -------------------------------

save = true                  # if true, saves output plots and run parameters
show = true                  # if true, shows plots

read_IR = true               # if true: read input IR pulse from file; if false: use Gaussian approximation 
read_ρ  = false              # if true: read gas density profile from file; if false: use pressure gradient approximation 
read_UV = false              # if true: overlay measured UV output on simulated results   

show_IR = false              # if true and "read_IR" is true: overlay measured input pulse on plots

# ----------------- INPUT HANDLING -------------------------------

# NOTE: see "./file_prepare.py" for relevant file content specifications

in_dir    = "input"                          # directory of input files

file_IR   = "IRpulse.dat"                    # name of IR input pulse file 
path_IR   = joinpath(in_dir, file_IR)        # sys. path to IR input pulse file 

file_ρ    = "dens.dat"                       # name of density profile data file 
path_ρ    = joinpath(in_dir, file_ρ)         # sys. path to density profile data file 

file_UV   = "UVpulse.dat"                    # name of IR input pulse file 
path_UV   = joinpath(in_dir, file_UV)        # sys. path to UV output pulse file 

# ------------------ SET MEASURED PARAMETERS ------------------------

# NOTE: some parameters below are overwritten when reading input data
#       directly from file (see QUICK SETTINGS above)

gas = :Ar           # gas
pres = 0.1          # central gas pressure [bar]
p_ed = 1e-3         # edge gas pressure [bar]
p_const = false     # if true: set constant pressure profile P==(pres,pres,pres) ; if false: set simple gradient: P==(p_ed, pres, p_ed)
τ = 5e-15           # FWHM pulse duration [s] (only relevant when temporal beam profile is approximated as Gaussian)
λ0 = 800e-9         # central wavelength [m]
w0 = 65e-6          # beam waist [m]
ϕ = 0.0             # central envelope offset (CEO) phase [rad]                                    -> can this be extracted from data?
energy = 300e-6     # pulse energy [J]                                                             -> multiply by 1kHz (?) repetition rate for beam power
L = 2e-3            # propagation distance (cell length) [m]

λ_lims = (200e-9, 1000e-9)      # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (200e-9, 360e-9)    # wavelength limits of UV region of interest [m,m]

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

R = 250e-6            # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  ->> Josina's version: 1e-3 
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

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),      # set nonlinear response of the gas: Kerr effect & plasma formation 
            Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))
            #Nonlinear.Kerr_env(PhysData.γ3_gas(gas)))

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

λhighidx    = argmin(abs.(λ .- λ_rangeUV[1]))         # index X such that λ[X]=λ_maxUV 
λlowidx     = argmin(abs.(λ .-λ_rangeUV[2]))          # index X such that λ[X]=λ_minUV

# * * * EXTRACT FIELD AMPLITUDES:
Eout = output.data["Eω"]                         # Hankel-transformed amplitude in frequency domain  
Erout = (q \ Eout)                               # Real-space amplitude in frequency domain at r ≠ 0   ->> "\" represents inverse Hankel transform 
Er0 = dropdims(Hankel.onaxis(Eout, q), dims=2)   # Real-space amplitude in frequency domain at r=0 

Etout = FFTW.irfft(Erout, length(t), 1)     # time-domain real field amplitude at r≠0
Et0 = FFTW.irfft(Er0, length(t),1)          # time-domain real field amplitude at r=0
Et = Maths.hilbert(Etout)                        # complex (analytic) representation of time-domain field amplitude 

# * * * CONVERT TO INTENSITIES:
Iωr = abs2.(Erout)                               #  intensity at r≠0 in frequency domain  [arbitrary units]
Iω0 = abs2.(Er0)                                 #  intensity at r=0 in frequency domain  [arbitrary units]
It0 = abs2.(Et0)                                 #  intensity at r=0 in time domain [arbitrary units]
It = PhysData.c * PhysData.ε_0/2 * abs2.(Et)     #  intensity in time domain in "physical units" [W/m²] 

# * * * FILTER FOR UV FIELD (r≠0):
filter=FFTW.rfft(Etout, 1)        # set up filter array 
filter[1:ωlowUVidx,:,:].=0;       # filters out ω < ω_min
filter[ωhighUVidx:end,:,:].=0;    # filters out ω > ω_max 

Etout_UV=FFTW.irfft(filter, length(t), 1)    # time-domain real field amplitude of UV pulse 

# * * * FILTER FOR UV FIELD (r=0)
filter_onaxis = FFTW.rfft(Et0, 1)          # set up filter array
filter_onaxis[1:ωlowUVidx,:].=0;           # filters out ω < ω_min
filter_onaxis[ωhighUVidx:end,:].=0;        # filters out ω > ω_max  

Et0_UV =FFTW.irfft(filter_onaxis, length(t), 1)    # time-domain real field amplitude of UV pulse at r=0
It0_UV = abs2.(Et0_UV)                           # intensity of on-axis UV pulse at 

# * * * EXTRACT INTENSITY ENVELOPES 
It0_envelope = abs2.(Maths.hilbert(Et0))          # envelope modulating It0
It0_UV_envelope = abs2.(Maths.hilbert(Et0_UV))    # envelope modulating It0_UV 

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

# * * * PROCESS MEASURED INPUT DATA 
if (read_IR & show_IR) == true
    
    # ADD CODE TO CONVERT TO FREQUENCY DOMAIN
end

if read_UV == true 
    λ_in    = readdlm(path_UV,' ', Float64, '\n')[:,1]                      # read in UV wavelength data
    I_in_UV = readdlm(path_UV,' ', Float64, '\n')[:,2]                      # read in UV intensity data 
    
    # ADD CODE TO CONVERT TO TIME DOMAIN

end    

# ----------------- OUTPUT HANDLING -------------------------

if isdir("output")== false   # create general directory to store output
    mkdir("output")
end    

if save == true              # create specific directory to store run output 
    out_path = joinpath("output","run_"*Dates.format(Dates.now(), "yyyy_mm_dd__HH_MM"))
    mkdir(out_path)
end 

# ----------------- PLOT RESULTS ----------------------------
import PyPlot:pygui, plt
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
    plt.savefig(joinpath(out_path,"off-axis_intensity.png"))
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
    plt.savefig(joinpath(out_path,"on-axis_intensity.png"))
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
    plt.savefig(joinpath(out_path,"density_and_susceptibility.png"))
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

if (read_IR & show_IR) == true 
    # ADD CODE TO OVERLAY INPUT IR SPECTRUM 
end    

plt.legend()

if save==true
    plt.savefig(joinpath(out_path,"full_on-axis_spectrum.png"))
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

plt.legend()

if save==true
    plt.savefig(joinpath(out_path,"UV_on-axis_spectrum.png"))
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

plt.legend()


if save==true
    plt.savefig(joinpath(out_path,"on-axis_spectrum_log.png"))
end 

#+++++ PLOT 7:  pulse energies and efficiency 
plt.figure(figsize=[7.04, 5.28])
plt.suptitle("THG conversion efficiency: η="*string(round(η_THG, digits=5) *100)*"%")
plt.subplots_adjust(hspace=0.5)

plt.subplot(2,1,1)
plt.plot(zout.*1e3, tot_pulse_en.*1e6, label="ΔE=-$(round(Int64, tot_pulse_en[1]*1e6-tot_pulse_en[end]*1e6))μJ", color="red")
plt.xlabel("z (mm)")
plt.ylabel("E (μJ)")
plt.title("Total pulse energy")
plt.legend()

plt.subplot(2,1,2)
plt.plot(zout.*1e3, UV_pulse_en.*1e9, label="ΔE=+$(round(Int64, UV_pulse_en[end]*1e9-UV_pulse_en[1]*1e9))nJ", color="red")
plt.xlabel("z (mm)")
plt.ylabel("E (nJ)")
plt.title("UV pulse energy")
plt.legend()

if save==true
    plt.savefig(joinpath(out_path,"pulse_energies.png"))
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
    plt.savefig(joinpath(out_path,"on-axis_frequency_evolution.png"))
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
    plt.plot(t_in*1e15,maximum(It0_envelope[:,1]).*Maths.normbymax(I_in), color="grey")
end

plt.legend(loc="upper right")

if save==true
    plt.savefig(joinpath(out_path,"time_domain_input.png"))
end

#+++++ PLOT 10: time-domain plot of UV output pulse 
plt.figure(figsize=[7.04, 5.28]) 
plt.title("Time-domain representation of on-axis UV output pulse")
plt.xlabel("t (fs)")
plt.ylabel("I(t; r=0, z=L) (arb. units)")
plt.plot(t*1e15, It0_UV[:,end] , color="red", label="FWHM pulse duration: τ="*string(round(τ_UV*1e15, digits=1) )*"fs")
plt.plot(t*1e15,It0_UV_envelope[:,end], color="black", linestyle="--")

if read_UV==true  # overlay measured UV output pulse 
    # ADD CODE TO OVERLAY UV PULSE IN TIME DOMAIN
end

plt.legend(loc="upper right")

if save==true
    plt.savefig(joinpath(out_path,"time_domain_UV_output.png"))
end

if show==true
    plt.show()
end    

# ----------------- WRITE PARAMS TO FILE ----------------------------

if save == true 
    open(joinpath(out_path,"params.txt"), "w") do file
        write(file, "gas     = :"*string(gas)*"\n")
        write(file, "pres    = "*string(pres)*"\n")
        write(file, "p_ed    = "*string(p_ed)*"\n")
        write(file, "p_const = "*string(p_const)*"\n")
        write(file, "τ       = "*string(τ)*"\n")
        write(file, "λ0      = "*string(λ0)*"\n")
        write(file, "w0      = "*string(w0)*"\n")
        write(file, "ϕ       = "*string(ϕ)*"\n")
        write(file, "energy  = "*string(energy)*"\n")
        write(file, "L       = "*string(L)*"\n")
        
        if read_IR == true
            write(file, "\n")
            write(file, "IR input beam read from: "*string(path_IR)*"\n")
        end

        if read_ρ == true
            write(file, "\n")
            write(file, "Gas density distribution read from: "*string(path_ρ)*"\n")
        end

        if read_UV == true
            write(file, "\n")
            write(file, "UV output beam read from: "*string(path_UV)*"\n")
        end
    end
end