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

provided option to base input pulse on measured data ; provided 
option to feed in custom gas density profile (untested!) ; 

		
	- compare input and output with measured data
 
		=> MUST BE PUT INTO IF STATEMENTS, IF NOT CRASHES PROGRAM WHEN TRYING TO PLOT!!!
			-> add extra if statement to manually turn off overlay?

	
		+ NIR spectrum scatter plot at input compared to measured data
			-> from Speck.dat  [check if FFT of Ek.dat gives good results]
		
		+ measured NIR pulse scatter plot at input compared to simulated NIR pulse at input 
			-> from Ek.dat 

		+ UV spectrum scatter plot at output compared to measured data
			-> from data folder [spectra ]

		+ measured UV pulse scatter plot at output compared to simulated UV pulse at output 
			-> from spectrum (provided FFT is good enough...)

	
	TAKE A STEP BACK AND REPRODUCE JOSINA'S ORIGINAL SIMULATIONS [SOME DISCREPANCIES!!!]

	
	PREPARE PRESENTATION FOR JOSINA


Stage 3: TO DO
-----------------------------

	- make simulation more realistic: 

		+ set up non-Gaussian spatial profile 

		+ input measured CEO phase value ?

	- produce gas density simulations	

	- vary pressure, beam intensity [when changing beam energy, also change beam shape (file)??], etc. to study effect on THG efficiency [Luna parameter scans...] (Ar & Ne)
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)
			-->> conversion efficiency and pulse energy!

