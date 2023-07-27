DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; 
	  carry out parameter scans to find optimal input beam 
          conditions for desired UV output beam;   
		 


Stage 1: DONE [Week 1-2 ]
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

Stage 2: IN PROGRESS [Week 2 - ?]
----------------------------------

	- use measured data as input (replacing "Fields.GaussGaussField") [DONE!]


	DOES INCREASING ENERGY HAVE DESIRED EFFECT (or do you have to also change input beam data accordingly?)

		
	- compare input and output with measured data 

		CHANGE COLOUR SCHEME ACCORDINGLY ? [RED, GREY, BLACK ?]

		+ NIR spectrum scatter plot at input compared to measured data
			-> from Speck.dat  
		+ measured NIR pulse scatter plot at input compared to simulated NIR pulse at input 
			-> from Ek.dat 

		+ UV spectrum scatter plot at output compared to measured data
			-> from data folder [spectra]
		+ UV pulse energy scatter plot compared to measured data
			-> from data folder [pressure scans]

	- DENSITY FUNCTION:
		literally just provide dens(z) [set up option to read from data file; cols: z[m], rho [1/m^2]] ]  

			γ = sellmeier_gas(gas)
			coren(ω; z) = sqrt(1 + γ(wlfreq(ω)*1e6)*dens(z))

	TAKE A STEP BACK AND REPRODUCE JOSINA'S ORIGINAL SIMULATIONS [SOME DISCREPANCIES!!!]

	
	PREPARE PRESENTATION FOR JOSINA


Stage 3: TO DO
-----------------------------

	- make simulation more realistic: 

		+ set up non-Gaussian spatial profile 

		+ input measured CEO phase value 

	- vary pressure, beam intensity, etc. to study effect on THG efficiency [Luna parameter scans...] (Ar & Ne)
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

