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

added different options for nonlinear model ; added visualisation
of pulse self-steepening; added option to overlay input pulse 
spectrum; 


- improve density distribution model
	+ wait for COMSOL data/access 
	+ change length/geometry of cell?

- adjust code
	+ improve structure for parameter scans (add loop option? read from param files for pressure scans ?)
	+ get access to Maxwell
		
- improve nonlinear effects
	+ TEST DIFFERENT OPTIONS 
		WHICH ONE IS BEST????!!!!    => JUST GO WITH ff FOR NOW!
	+ read papers 

- carry out parameter scans and see if filamentation effects (saturation, 
  second spectral peak, beam deformation) set in 
	
	+ vary pressure
	+ vary gas 
	+ vary beam energy 
	+ vary CEO phase 
	+ turn ionisation on/off

	=> find conditions to maximise conversion efficiency (and 
	   minimise pulse duration?)





