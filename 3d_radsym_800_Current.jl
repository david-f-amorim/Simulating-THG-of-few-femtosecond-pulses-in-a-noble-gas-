using Luna
import Luna.PhysData: wlfreq
import FFTW
import Luna: Hankel
import NumericalIntegration: integrate, SimpsonEven

# ------------------ SET MEASURED PARAMETERS ------------------------

gas = :Ar           # gas
pres = 0.4          # central gas pressure [bar]
p_ed = 1e-3         # edge gas pressure [bar]
τ = 5e-15           # FWHM pulse duration [s]
λ0 = 800e-9         # central wavelength [m]
w0 = 65e-6          # beam waist [m]
energy = 150e-6     # pulse energy [J]
L = 2e-3            # propagation distance (cell length) [m]

λ_lims = (200e-9, 1000e-9)      # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (200e-9, 360e-9)    # wavelength limits of UV region of interest [m,m]

# ----------------- SET PRESSURE PROFILE ---------------------------

Z= (0, L/2, L)        # define points for pressure gradient
P= (p_ed, pres, p_ed) # define values for pressure gradient [constant pressure: P==(pres,pres,pres); simple gradient: P==(p_ed, pres, p_ed)]

(coren,dens)=Capillary.gradient(gas,Z,P)   # gives n(ω; z) and ρ(z) from pressure profile ->> IS THIS UPDATED TO ACCOUNT FOR IONISATION?

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

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),    # set nonlinear response of the gas: Kerr effect & plasma formation
            Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))       # ->> "grid.to" appears twice; BUG OR FEATURE ??

# ----------------- SET INPUT FIELD ----------------------------                                       

inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ,             # models input beam as Gaussian (TO BE REPLACED WITH "DataField" BASED ON MEASURED SPECTRUM)
            energy=energy, w0=w0, propz=-0.001)

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
ω0idx       = argmin(abs.(grid.ω .- 2π*PhysData.c/λ0))                   # index X such that ω[X]=ω0 (fundamental)
ωTHidx      = argmin(abs.(grid.ω .- 2π*PhysData.c/(λ0/3)))               # index X such that ω[X]=3ω0 (TH)
ωhighUVidx    = argmin(abs.(grid.ω .- 2π*PhysData.c/λ_rangeUV[1]))         # index X such that ω[X]=ω_maxUV 
ωlowUVidx     = argmin(abs.(grid.ω .- 2π*PhysData.c/λ_rangeUV[2]))         # index X such that ω[X]=ω_minUV

λ0idx       = argmin(abs.(λ .- λ0))                   # index X such that λ[X]=λ0 (fundamental)
λTHidx      = argmin(abs.(λ .- (λ0/3)))               # index X such that λ[X]=λ/3 (TH))
λhighidx    = argmin(abs.(λ .- λ_rangeUV[1]))         # index X such that λ[X]=λ_maxUV 
λlowidx     = argmin(abs.(λ .-λ_rangeUV[2]))          # index X such that λ[X]=λ_minUV

f0 = f[ω0idx]     # linear frequency of fundamental
fTH = f[ωTHidx]   # linear frequency of third harmonic 

# * * * EXTRACT FIELD AMPLITUDES:
Eout = output.data["Eω"]                         # Hankel-transformed amplitude in frequency domain  
Erout = (q \ Eout)                               # Real-space amplitude in frequency domain at r ≠ 0   ->> "\" represents inverse Hankel transform 
Er0 = dropdims(Hankel.onaxis(Eout, q), dims=2)   # Real-space amplitude in frequency domain at r=0 

Etout = FFTW.irfft(Erout, length(grid.t), 1)     # time-domain real field amplitude 
Et = Maths.hilbert(Etout)                        # complex (analytic) representation of time-domain field amplitude 

# * * * CONVERT TO INTENSITIES:
Iωr = abs2.(Erout)                               #  intensity at r≠0 in frequency domain  [arbitrary units]
Iω0 = abs2.(Er0)                                 #  intensity at r=0 in frequency domain  [arbitrary units]
It = PhysData.c * PhysData.ε_0/2 * abs2.(Et)     #  intensity in time domain in "physical units" [W/m²] 

# * * * FILTER FOR UV FIELD:
filter=FFTW.rfft(Etout, 1)        # set up filter array 
filter[1:ωlowUVidx,:,:].=0;       # filters out ω < ω_min
filter[ωhighUVidx:end,:,:].=0;    # filters out ω > ω_max 

Etout_UV=FFTW.irfft(filter, length(grid.t), 1)    # time-domain real field amplitude of UV pulse 

# * * * CALCULATE PULSE ENERGIES:
tot_pulse_en = zeros(length(zout))                   # set up array for total pulse energy 
UV_pulse_en  = zeros(length(zout))                   # set up array for UV pulse energy 

for i = 1:size(Etout, 3)
    tot_pulse_en[i] = energyfun(Etout[:, :, i])         # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
    UV_pulse_en[i]  = energyfun(Etout_UV[:, :, i])      # fill array: tot_pulse_en[i] is the total pulse energy at some "z" value
end     

η_THG = UV_pulse_en[end]/tot_pulse_en[1]            # THG efficiency: initial total pulse energy / final UV pulse energy 

# ----------------- PLOT RESULTS ----------------------------
import PyPlot:pygui, plt
close("all")
pygui(true)

#+++++ PLOT 1:  fundamental and third harmonic intensities as functions of z and r≠0
plt.figure()
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
    
#+++++ PLOT 2:  fundamental and third harmonic intensities as functions of z at r=0
plt.figure()
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

#+++++ PLOT 3: gas number density and effective susceptibility along the cell 
plt.figure()
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

#+++++ PLOT 4:  linear on-axis spectrum I(λ) at z=0 and z=L 
plt.figure()
plt.suptitle("Linear on-axis beam spectrum (total & UV-only)")
plt.subplots_adjust(hspace=0.3)

plt.subplot(2,1,1)
plt.plot(λ[2:end]*1e9, Iω0[2:end,1], label="z=0mm", color="grey")
plt.plot(λ[2:end]*1e9, Iω0[2:end,end], label="z=$(L*1e3)mm", color="red")
plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ) (arb. units)")
plt.legend()

plt.subplot(2,1,2)
plt.plot(λ[λlowidx:λhighidx]*1e9, Iω0[λlowidx:λhighidx,1], label="z=0mm", color="grey")
plt.plot(λ[λlowidx:λhighidx]*1e9, Iω0[λlowidx:λhighidx,end], label="z=$(L*1e3)mm", color="red")
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ) (arb. units)")
plt.legend()

#+++++ PLOT 5:  log. on-axis spectrum I(λ) at z=0 and z=L 
Iω0log = log10.(Iω0)

plt.figure()
plt.suptitle("Log. on-axis beam spectrum (total & zoomed)")
plt.subplots_adjust(hspace=0.3)

plt.subplot(2,1,1)
plt.plot(λ[2:end]*1e9, Iω0log[2:end,1], label="z=0mm",  color="grey")
plt.plot(λ[2:end]*1e9, Iω0log[2:end,end], label="z=$(L*1e3)mm", color="red")
plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ) (arb. units)")
plt.legend()

plt.subplot(2,1,2)
plt.plot(λ[2:end]*1e9, Iω0log[2:end,1], label="z=0mm", color="grey")
plt.plot(λ[2:end]*1e9, Iω0log[2:end,end], label="z=$(L*1e3)mm", color="red")
plt.xlim(λ_lims[1]*1e9, λ_lims[2]*1e9)
plt.ylim(bottom=10)
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ) (arb. units)")
plt.legend()

#+++++ PLOT 6:  pulse energies and efficiency 
plt.figure()
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

#+++++ PLOT 7: on-axis frequency evolution 
plt.figure()
plt.suptitle("On-axis frequency evolution")
plt.pcolormesh(zout*1e3, f*1e-15, log10.(Maths.normbymax(Iω0)))    # originally: logarithmic 
plt.colorbar()
plt.clim(0, -6)  # cut away tails of distr. 
plt.xlabel("z (mm)")
plt.ylabel("f (PHz)")
plt.title("log. I(r=0, ω), normed")

# !!! fix this



plt.show()