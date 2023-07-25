GOAL: improve THG simulation
-----------------------------

Stage 1: 

	- fix existing code for better comparison with Josina's:

		+ fix efficiency display rounding problem
		+ add normed LOG spectrum / norm existing log spectrum ?
		+ fix frequency evolution: Josina normalised before taking log?? add legend/arrows for f0 and fTH

		-> try to reproduce Josina's results!!
	
		
	- produce additional output visualisations: 
		

		+ phase information ?? 
		+ time-domain representation (including TL pulse length & time-bandwidth product??["Processing.jl"]) ??


Stage 2:

	- set up "Fields.DataField" to represent input pulse (replacing "Fields.GaussGaussField")
		
		+ extract spectral phase from FROG data  (read through Jupyter Notebook)
		+ combine spectral phase with frequency and spectral power density
                  into a data file 

	- compare input and output with measured data 

		+ NIR spectrum at input compared to measured data
		+ UV spectrum at output compared to measured data
		+ ....


Stage 3:

	- vary pressure, beam intensity, etc. to study effect on THG efficiency
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

	- set up own density function