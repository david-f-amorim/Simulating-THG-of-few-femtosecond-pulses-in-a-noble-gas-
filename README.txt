GOAL: improve THG simulation
-----------------------------

Step 1: 

	- reproduce Josina's existing simulation results (familiarisation with code)
	- tidy up code where necessary

		+ ensure origin of all variable values are clear
		+ fix input variables for "Nonlinear.PlasmaCumtrapz" in "responses" ?
		
	- produce suitable output visualisation for easy comparison with measured data

		+ refractive index as function of propagation distance: I(z)
		+ gas density as function of propagation distance: rho(z)
		+ UV energy [as func of z and total] and conversion efficiency
		+ [log and non log] spectrum at z=0 and z=L: I(f) or I(lam)
		+ NIR spectrum at input and output [log and non-log; overlay of measured data?]
		+ UV spectrum at output [log and non-log; overlay of measured data?]
		+ beam profile I(r,z) for both NIR and UV 
		+ radial beam profile at z=0 for NIR as well as at z=L for NIR & UV 
		
		+ phase information ?? time-domain representation ?? pulse length ?? time-bandwidth product ?? ["Processing.jl"]

		+ INTERPRET CURRENT PLOTS IN THE CODE!!!
		  

Step 2:

	- set up "Fields.DataField" to represent input pulse (replacing "Fields.GaussGaussField")
		
		+ extract spectral phase from FROG data  (read through Jupyter Notebook)
		+ combine spectral phase with frequency and spectral power density
                  into a data file 