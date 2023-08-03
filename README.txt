DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; carry out parameter scans to find optimal 
	  input beam conditions for desired UV output beam;   
		 

Stage 1: DONE [Weeks 1-2]
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

Stage 2: DONE [Week 2]
----------------------------------

provided option to base input pulse on measured data ; provided 
option to feed in custom gas density profile (untested!) ; 
provided option to overlay measured UV spectrum 
	

Stage 3: IN PROGRESS [Weeks 3 -]
--------------------------------  

added different options for nonlinear model (use "ff" for now) ; 
added visualisation of pulse self-steepening; added option to overlay
input pulse spectrum; changed cell length to 3mm ; implemented option 
to run pressure scans ; 


- adjust code
	+ write code to process pressure scan results:
		ONLY PLOT UP TO 10-15 PRESSURES FOR SPECTRA!

		THINK ABOUT HOW TO BEST COMPARE OUTPUT FROM DIFFERENT SCANS!
			[ENERGY/EFFICIENCY]

		THINK ABOUT HOW TO COMPARE SPECTRA TO MEASURED DATA 

	+ get access to Maxwell
		
- carry out parameter scans and see if filamentation effects (saturation, 
  second spectral peak, beam deformation) set in 
	
	=> first generate lots of scan data for general analysis, 
	   then look at individual cases in more detail

	WHICH SCAN TO RUN NEXT???!!!
		Ar, f, 0.0, ion at different energies?	
	
	+ vary pressure 
	+ vary beam energy 
	+ vary CEO phase

	+ vary gas 
	+ turn ionisation on/off
	+ vary kerr 

	=> find conditions to maximise conversion efficiency (and 
	   minimise pulse duration?)

- do self-steepening investigations -> leading or trailing edge ??

- improve density distribution model
	+ wait for COMSOL data/access 

- improve nonlinear effects
	+ test different options => JUST GO WITH ff FOR NOW!
	+ read papers 


