using Luna
import Luna.PhysData: wlfreq
import FFTW
import Luna: Hankel
import NumericalIntegration: integrate, SimpsonEven


gas = :Ar           # gas
pres = 0.4          # pressure (bar)
low = 1e-3          # pressure on the outside of the pressure profile
τ = 5e-15           # pulse duration ?, full width half maximum
λ0 = 800*1e-9       # fundamental central wavelength
lam=λ0              # ?? --> not needed?

λlow = 200e-9       # low wavelength limit of UV region of interest (for plotting??, FFT???)
λhigh = 360e-9      # high wavelength limit of UV region of interest (for plotting??, FFT???)

w0 = 65e-6          # waist, definition: w(z) is the radius at which the field amplitudes fall to 1/e of their axial values (i.e., where the intensity values fall to 1/e2 of their axial values)
energy = 150e-6    # pulse energy
L = 0.002           # total distance to propagate

Z=(0, L/2, L) #defining points for presure gradient   
P=(pres, pres, pres) # values for pressure gradient
(coren,dens)=Capillary.gradient(gas,Z,P);

R = 1e-3            # aperture radius for Hankel transform    ## WHY ???
N = 1024            # sample size for Hankel tansform         ## WHY ???

grid = Grid.RealGrid(L, 800e-9, (200e-9, 1000e-9), 0.05e-12) 
q = Hankel.QDHT(R, N, dim=2)

energyfun, energyfun_ω = Fields.energyfuncs(grid, q)

densityfun = let dens0=PhysData.density(gas, pres)   # currently not used ?
    z -> dens0
end

ionpot = PhysData.ionisation_potential(gas)
ionrate = Ionisation.ionrate_fun!_PPTcached(gas, λ0)

responses = (Nonlinear.Kerr_field(PhysData.γ3_gas(gas)),Nonlinear.PlasmaCumtrapz(grid.to, grid.to, ionrate, ionpot)) # !!!!

#linop = LinearOps.make_const_linop(grid, q, PhysData.ref_index_fun(gas, pres))
linop = LinearOps.make_linop(grid, q, coren)
#normfun = NonlinearRHS.const_norm_radial(grid, q, PhysData.ref_index_fun(gas, pres))
normfun = NonlinearRHS.norm_radial(grid, q, coren)
inputs = Fields.GaussGaussField(λ0=λ0, τfwhm=τ, energy=energy, w0=w0, propz=-L/2)           # !!!!!!!!!!

#Eω, transform, FT = Luna.setup(grid, q, densityfun, normfun, responses, inputs)
Eω, transform, FT = Luna.setup(grid, q, dens, normfun, responses, inputs)

# statsfun = Stats.collect_stats(grid, Eω, Stats.ω0(grid))
output = Output.MemoryOutput(0, grid.zmax, 201)
Luna.run(Eω, grid, linop, transform, FT, output)



#-------------------------------------------------------------------------------------------


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
energy = zeros(length(zout))
for ii = 1:size(Etout, 3)
    energy[ii] = energyfun(Etout[:, :, ii]);
    Ir[:, ii] = integrate(grid.ω, Iωr[:, :, ii], SimpsonEven());
end

ω0idx = argmin(abs.(grid.ω .- 2π*PhysData.c/λ0))
ω133idx = argmin(abs.(grid.ω .- 2π*PhysData.c/133e-9))
ωhighidx= argmin(abs.(grid.ω .- 2π*PhysData.c/λlow))
ωlowidx= argmin(abs.(grid.ω .- 2π*PhysData.c/λhigh))

zr = π*w0^2/λ0
points = L/2 .+ [-15, 3, 21].*zr
idcs = [argmin(abs.(zout .- point)) for point in points]

Epoints = [Hankel.symmetric(Et[:, :, idxi], q) for idxi in idcs]
rsym = Hankel.Rsymmetric(q);

tofilter=FFTW.rfft(Etout, 1)
tofilter[1:ωlowidx,:,:].=0;
tofilter[ωhighidx:end,:,:].=0;
back=FFTW.irfft(tofilter, length(grid.t), 1)


energy2 = zeros(length(zout));
for ii = 1:length(zout)
    energy2[ii] = energyfun(back[:, :, ii]);
#   Iyx[:, :, ii] = (grid.ω[2]-grid.ω[1]) .* sum(Iωyx[:, :, :, ii], dims=1);
end

efficiency=energy2[end]/energy[1]

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

plt.figure()
plt.subplot(2,1,1)
plt.pcolormesh(zout, q.r, Iωr[ω0idx, :, :])
plt.colorbar()
plt.ylim(0, 0.0002)
plt.ylabel("r (m)")
plt.title("I(r, ω=400nm)")

plt.subplot(2,1,2)
plt.pcolormesh(zout, q.r, Iωr[ω133idx, :, :])
plt.colorbar()
plt.ylim(0, 0.0002)
plt.xlabel("z (m)")
plt.ylabel("r (m)")
plt.title("I(r, ω=133nm)")

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

# plt.figure()
# plt.pcolormesh(zout*1e2, grid.ω*1e-15/2π, Iω0log)
# plt.colorbar()
# plt.clim(0, -6)
# plt.xlabel("z (m)")
# plt.ylabel("f (PHz)")
# plt.title("I(r=0, ω)")

# plt.figure()
# plt.plot(zout.*1e2, energy.*1e6)
# plt.xlabel("Distance [cm]")
# plt.ylabel("Energy [μJ]")

# jw = Plotting.cmap_white("jet"; n=10)
# fig = plt.figure()
# fig.set_size_inches(12, 4)
# for ii in 1:3
#     plt.subplot(1, 3, ii)
#     plt.pcolormesh(grid.t*1e15, rsym*1e3, abs2.(Epoints[ii]'), cmap=jw)
#     plt.xlim(-42, 42)
#     plt.ylim(-1.8, 1.8)
# end

plt.figure()
plt.plot(zout.*1e2, energy2.*1e6)
plt.xlabel("Distance (cm)")
plt.ylabel("VUV pulse energy (μJ)")
#plt.title("Wavelength = "*string(λ0)*"; Efficiency = "*string(efficiency))

# save_results_to ="Z:\\user\\Nora Schmitt\\Numerics\\julia\\scan-images\\"
# plt.savefig(save_results_to*string(λ0)*".png", dpi = 300)

profile_133=Iωr[ω133idx,:,:]
for i=1:length(zout)
profile_133[:,i]=Iωr[ω133idx,:,length(zout)]
end

fig = plt.figure()
ax = fig.add_subplot(111, projection="polar")
θ = LinRange(0., 2π, length(zout))
pc = plt.pcolormesh(θ,q.r,profile_133)
cbar = plt.colorbar(pc)
cbar.set_label("Intensity")
ax[:grid](true)
ax[:set_yticks]([0, 0.00005])
ax[:set_yticklabels]([" ","50"])

plt.show()
