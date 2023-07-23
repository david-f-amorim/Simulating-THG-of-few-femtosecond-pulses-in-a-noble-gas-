using Luna
import Luna.PhysData: wlfreq
import FFTW
import Luna: Hankel
import NumericalIntegration: integrate, SimpsonEven

# ------------------ SET MEASURED PARAMETERS ------------------------

gas = :Ar           # gas
pres = 0.4          # central gas pressure [bar]
pres_edge = 1e-3    # edge gas pressure [bar]
τ = 5e-15           # FWHM pulse duration [s]
λ0 = 800e-9         # central wavelength [m]
w0 = 65e-6          # beam waist [m]
energy = 150e-6     # pulse energy [J]
L = 2e-3            # propagation distance (cell length) [m]

λ_lims = (200e-9, 1000e-9)      # wavelength limits of overall frequency window (both NIR and UV) [m,m]
λ_rangeUV = (200e-9, 360e-9)    # wavelength limits of UV region of interest [m,m]

# ----------------- SET PRESSURE PROFILE ---------------------------
#                                         Note: model could be improved?!

Z= (0, L/2, L)        # define points for pressure gradient
P= (pres, pres, pres) # define values for pressure gradient [constant pressure: P==(pres,pres,pres); simple gradient: P==(pres_edge, pres, pred_edge)]

(coren,dens)=Capillary.gradient(gas,Z,P)   # gives n(z) and ρ(z) from pressure profile 

# ----------------- SET SIMULATION GRID ----------------------------

R = 250e-6            # aperture radius for Hankel transform (assume field ≈0 for r>R ) [m]  ->> should be 250e-6 m ? ; Josina's version: 1e-3 ; AFFECTS RUNTIME?? [should be along r not z?!]
N = 1024              # sample size for Hankel tansform grid     [1]                         ->> WHY THIS VALUE ?? WOULD CHOOSING N≫ IMPROVE THINGS?
trange = 0.05e-12     # total extent of time window required [s]                             ->> Josina's version: 0.05e-12 ; SHOULD BE INCREASED TO trange = L/c ∼ 1e-11 ???


grid = Grid.RealGrid(L, λ0, λ_lims ,trange)   # set up time grid
q = Hankel.QDHT(R, N, dim=2)                  # set up discrete Hankel transform matrix      ->> dim represents axis of transform; WHICH ONE?? k or z???                    

energyfun, _ = Fields.energyfuncs(grid, q)    #                                               ->> ??????


# ----------------- SET NONLINEAR EFFECTS ----------------------------

ionpot = PhysData.ionisation_potential(gas)          
ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)

linop = LinearOps.make_linop(grid, q, coren)
normfun = NonlinearRHS.norm_radial(grid, q, coren)

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),
            Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot))  # ->> "grid.to" appears twice FIX THIS ?!

# ----------------- SET INPUT FIELD ----------------------------
#                                       replace with DataField!

inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ, energy=energy, w0=w0, propz=-0.001)

# ----------------- RUN SIMULATION ----------------------------

Eω, transform, FT = Luna.setup(grid, q, dens, normfun, responses, inputs)
output = Output.MemoryOutput(0, grid.zmax, 201)    # ???!!!!
Luna.run(Eω, grid, linop, transform, FT, output)            


# -----------------------------------------------------------------------------------------------------
#     PLOTTING AND PROCESSING RESULTS
# -----------------------------------------------------------------------------------------------------

ω = grid.ω
t = grid.t

zout = output.data["z"]
Eout = output.data["Eω"]

Erout = (q \ Eout)
Iωr = abs2.(Erout)
# Iω0 = Iωr[:, 1, :]
Er0 = dropdims(Hankel.onaxis(Eout, q), dims=2);
Iω0 = abs2.(Er0);
Iω0log = log10.(Maths.normbymax(Iω0));
Etout = FFTW.irfft(Erout, length(grid.t), 1)

Ilog = log10.(Maths.normbymax(abs2.(Eout)))

It = PhysData.c * PhysData.ε_0/2 * abs2.(Maths.hilbert(Etout));
Itlog = log10.(Maths.normbymax(It))

Ir = zeros(Float64, (length(q.r), length(zout)))

Et = Maths.hilbert(Etout)
energy = zeros(length(zout))                    # IR (?) (or total?) pulse energy
for ii = 1:size(Etout, 3)
    energy[ii] = energyfun(Etout[:, :, ii]);
    Ir[:, ii] = integrate(grid.ω, Iωr[:, :, ii], SimpsonEven());
end

ω0idx = argmin(abs.(grid.ω .- 2π*PhysData.c/λ0))
ω133idx = argmin(abs.(grid.ω .- 2π*PhysData.c/133e-9))
ωhighidx= argmin(abs.(grid.ω .- 2π*PhysData.c/λ_rangeUV[1]))
ωlowidx= argmin(abs.(grid.ω .- 2π*PhysData.c/λ_rangeUV[2]))

zr = π*w0^2/λ0
points = L/2 .+ [-15, 3, 21].*zr
idcs = [argmin(abs.(zout .- point)) for point in points]

Epoints = [Hankel.symmetric(Et[:, :, idxi], q) for idxi in idcs]
rsym = Hankel.Rsymmetric(q);

tofilter=FFTW.rfft(Etout, 1)
tofilter[1:ωlowidx,:,:].=0;
tofilter[ωhighidx:end,:,:].=0;
back=FFTW.irfft(tofilter, length(grid.t), 1)        # 256×1024×201 Array


energy2 = zeros(length(zout));                      # UV pulse energy
for ii = 1:length(zout)
    energy2[ii] = energyfun(back[:, :, ii]);
#   Iyx[:, :, ii] = (grid.ω[2]-grid.ω[1]) .* sum(Iωyx[:, :, :, ii], dims=1);
end

efficiency=energy2[end]/energy[1]

########################################################################################################################### plotting results

import PyPlot:pygui, plt
close("all")
pygui(true)
Iλ0 = Iωr[ω0idx, :, :]
#λ0 = 2π*PhysData.c/grid.ω[ω0idx]
w1 = w0*sqrt(1+(L/2*λ0/(π*w0^2))^2)
# # w1 = w0
# Iλ0_analytic = Maths.gauss.(q.r, w1/2)*(w0/w1)^2 # analytical solution (in paraxial approx)
# plt.figure()
# plt.plot(q.r*1e3, Maths.normbymax(Iλ0[:, end]))
# plt.plot(q.r*1e3, Maths.normbymax(Iλ0_analytic), "--")

## plot third harmonic and fundamental intensity along x and r
plt.figure()
plt.tight_layout()
plt.subplot(2,1,1)
plt.tight_layout()
plt.pcolormesh(zout*1e3, q.r*1e3, Iωr[ω0idx, :, :])
plt.colorbar()
#plt.ylim(0, 0.2)
plt.ylabel("r (mm)")
plt.title("I(r, ω=800nm)")

plt.subplot(2,1,2)
plt.tight_layout()
plt.pcolormesh(zout*1e3, q.r*1e3, Iωr[ω133idx, :, :])
plt.colorbar()
plt.ylim(0, 0.2)
plt.xlabel("z (mm)")
plt.ylabel("r (mm)")
plt.title("I(r, ω=266nm)")
plt.tight_layout()

# plt.figure()
# plt.pcolormesh(zout, q.r, It[length(grid.t)÷2, :, :])
# plt.colorbar()
# plt.ylim(0, 4)
# plt.xlabel("z (m)")
# plt.ylabel("r (m)")
# plt.title("I(r, t=0)")

# plt.figure()
# plt.pcolormesh(zout*1e2, q.r*1e3, Ir)
# # plt.pcolormesh(zout*1e2, q.r*1e3, log10.(Maths.normbymax(Ir)))
# plt.colorbar()
# plt.ylim(0, R*1e3)
# # plt.clim(-6, 0)
# plt.xlabel("z (m)")
# plt.ylabel("r (m)")
# plt.title("\$\\int I(r, \\omega) d\\omega\$")

## plot frequency propagation
plt.figure()
plt.pcolormesh(zout*1e2, grid.ω*1e-15/2π, Iω0log)
plt.colorbar()
plt.clim(0, -6)
plt.xlabel("z (m)")
plt.ylabel("f (PHz)")
plt.title("I(r=0, ω)")

## plot wavelength propagation
# plt.figure()
# plt.pcolormesh(zout*1e2, (2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0log[2:end,:])
# plt.colorbar()
# plt.clim(0, -6)
# plt.xlabel("z (m)")
# plt.ylabel("λ (nm)")
# plt.title("I(r=0, λ)")

## plot final spectrum, frequency
plt.figure()
plt.tight_layout()
# plt.plot(grid.ω*1e-15/2π, Iω0[:,1], label="z=0mm", color="grey")
plt.plot(grid.ω*1e-15/2π, Iω0[:,1]./maximum(Iω0[:,1]), label="z=0mm", color="grey")        # normalized
# plt.plot(grid.ω*1e-15/2π, Iω0[:,end], label="z=2mm")
plt.plot(grid.ω*1e-15/2π, Iω0[:,end]./maximum(Iω0[:,end]), label="z=2mm")                  # normalized
# plt.ylim(-6, 0)
plt.xlabel("f (PHz)")
plt.ylabel("I(r=0, f), norm.")
plt.legend()

## plot final spectrum, wavelength
plt.figure()
plt.tight_layout()
plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0[2:end,1], label="z=0mm", color="grey")
# plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0[2:end,1]./maximum(Iω0[2:end,1]), label="z=0mm", color="grey")     # normalized
plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0[2:end,end], label="z=2mm")
# plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0[2:end,end]./maximum(Iω0[2:end,end]), label="z=2mm")               # normalized
# plt.xlim(150, 1500)
plt.xlim(180, 1250)
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ)")
plt.legend()


## plot final spectrum, frequency, log
plt.figure()
plt.tight_layout()
# plt.plot(grid.ω*1e-15/2π, Iω0log[:,1], label="z=0mm", color="grey")
plt.plot(grid.ω*1e-15/2π, -1*Iω0log[:,1]./minimum(Iω0log[:,1]), label="z=0mm", color="grey")            # normalized
# plt.plot(grid.ω*1e-15/2π, Iω0log[:,end], label="z=2mm")
plt.plot(grid.ω*1e-15/2π, -1*Iω0log[:,end]./minimum(Iω0log[:,end]), label="z=2mm")                      # normalized
# plt.ylim(-6, 0)
plt.xlabel("f (PHz)")
plt.ylabel("I(r=0, f), norm.")
plt.legend()

## plot final spectrum, wavelength, log
plt.figure()
plt.tight_layout()
plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0log[2:end,1], label="z=0mm", color="grey")
# plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], -1*Iω0log[2:end,1]./minimum(Iω0log[2:end,1]), label="z=0mm", color="grey")     # normalized
plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], Iω0log[2:end,end], label="z=2mm")
# plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end], -1*Iω0log[2:end,end]./minimum(Iω0log[2:end,end]), label="z=2mm")               # normalized
# plt.xlim(150, 1500)
plt.xlim(180, 1250)
plt.ylim(-10, 0)
plt.xlabel("λ (nm)")
plt.ylabel("I(r=0, λ)")
plt.legend()

## plot refractive index
# plt.figure()
# plt.tight_layout()
# # plt.plot(grid.ω*1e-15/2π, coren.(grid.ω; z=0.001), label="refractive index of Argon at z=0.001", color="grey") #frequency axis
# # plt.xlabel("f (PHz)")
# # plt.ylabel("refractive index")
# plt.plot((2π*PhysData.c*ones(length(grid.ω))./grid.ω *1e9)[2:end],coren.(grid.ω; z=0.001)[2:end], label="refractive index of Argon at z=0.001", color="grey") #wavelength axis
# plt.xlabel("λ (nm)")
# plt.ylabel("refractive index")


## IR (?) (or total?) pulse energy along the propagation axis
# plt.figure()
# plt.plot(zout.*1e2, energy.*1e6)
# plt.xlabel("Distance [cm]")
# plt.ylabel("Energy [μJ]")

jw = Plotting.cmap_white("jet"; n=10)
fig = plt.figure()
fig.set_size_inches(12, 4)
for ii in 1:3
    plt.subplot(1, 3, ii)
    plt.pcolormesh(grid.t*1e15, rsym*1e3, abs2.(Epoints[ii]'), cmap=jw)
    plt.xlim(-42, 42)     # time (fs)
    plt.ylim(-1.8, 1.8)   # rsym (mm)
end

## UV pulse energy along the propagation axis
plt.figure()
plt.plot(zout.*1e3, energy2.*1e9, label="pulse energy")
plt.xlabel("Distance (mm)")
plt.ylabel("Pulse energy (nJ)")
#plt.vlines(1, 0, 800, color="red", label="focus position")
#plt.legend()
#plt.title("Wavelength = "*string(λ0)*"; Efficiency = "*string(efficiency))

# save_results_to ="Z:\\user\\Nora Schmitt\\Numerics\\julia\\scan-images\\"
# plt.savefig(save_results_to*string(λ0)*".png", dpi = 300)

## plot polar UV beam shape
# profile_133=Iωr[ω133idx,:,:]
# for i=1:length(zout)
# profile_133[:,i]=Iωr[ω133idx,:,length(zout)]
# end
# fig = plt.figure()
# ax = fig.add_subplot(111, projection="polar")
# θ = LinRange(0., 2π, length(zout))
# pc = plt.pcolormesh(θ,q.r,profile_133)
# cbar = plt.colorbar(pc)
# cbar.set_label("Intensity")
# ax[:grid](true)
# ax[:set_yticks]([0, 0.00005])
# ax[:set_yticklabels]([" ","50"])

#= COMMENTED OUT BC OF PYTHON ERROR; FIX LATER!

## Plot gas density
dens_x = zeros(0)
for x in zout
    append!(dens_x, dens(x))
end
plt.figure()
plt.tight_layout()
plt.plot(zout.*1e3, dens_x)
plt.xlabel("Distance (mm)")
plt.ylabel("density (kg/m3)")
plt.pcolormesh(zout, q.r, dens_x)
plt.colorbar()
plt.ylim(0, 0.0002)
plt.xlabel("z (m)")
plt.ylabel("r (m)")
plt.title("density")
 =#
plt.show()


## plot experimental data
# using Glob

# folderPath = readdir("A://data//2022//2022-01-06//Spectra//UV/Argon//Ar Pressure scan 150mW IR", join=true)

# data_dict = Dict{String, Array{Float64, 2}}()

# for fileName in folderPath
#     if last(fileName, 4) == ".txt"
#         lines = readlines(fileName)
#         index_BeginSpectralData = findfirst(x -> contains(x, ">Begin Spectral Data<"), lines)[1] + 1
#         # lines = [split(line, '\t') for line in lines]
#         data = [parse(Float64, x) for x in lines[index_BeginSpectralData:end]]
#         data[:, 2] ./= maximum(data[:, 2])  # normalizes each (!) spectrum
#         close(fileName)
#         pressure = split(split(fileName, '\\')[end], '_')[1]
#         data_dict[pressure] = data
#     end
# end

# h = 4.14e-15  # eV*s
# c = 2.99792458e17  # nm/s

# pressures_to_plot = [
# #    "0.1bar",
# #    "0.2bar",
# #    "0.3bar",
#     "0.4bar",
# #    "0.5bar",
# #    "0.6bar",
# #    "0.7bar",
# #    "0.8bar",
# #    "1.0bar",
# #    "1.2bar",
# #    "1.4bar",
# #    "1.5bar",
# #    "1.6bar",
# #    "2.0bar",
# #    "2.5bar",
# #    "3.0bar",
# #    "3.5bar",
# #    "4.0bar",
# #    "4.5bar",
# #    "5.0bar"
# ]

# #evenly_spaced_interval = range(1, stop=0, length=length(data_dict))
# #colors = plt.colormap("plasma")(evenly_spaced_interval)
# figure(figsize=(6, 6), dpi=80)

# for (i, line) in enumerate(sort(collect(data_dict), rev=false))
#     pressure = line[1]
#     if pressure in pressures_to_plot
#         println(line)
#         lambd = line[2][:, 1]
#         energy = h * c ./ lambd
# #        x = energy
#         x = lambd
#         plt.plot(x, line[2][:, 2], label=pressure, linestyle="-", color="blue")
# #        plt.fill_between(x, line[2][:, 2], label=pressure, linestyle="-", facecolor="blue", alpha=0.4)
#     end
# end

# tick_params(direction="in", top=true, right=true)
# xlabel("Wavelength (nm)")
# # xlabel("Photon energy (eV)")
# ylabel("Intensity (norm.)")
# tight_layout()


