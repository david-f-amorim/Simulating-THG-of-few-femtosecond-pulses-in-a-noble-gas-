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

- improve density distribution model
	+ wait for COMSOL data/access 
	+ qualitatively compare existing COMSOL 
	  plots with gradient approximation 
	+ change length/geometry of cell?

- adjust code
	+ include option to read in spectral IR data for comparison with spectrum
	  calculated from timing data 
	+ improve structure for parameter scans (add loop option?)
	+ get access to Maxwell?
	+ add plot temporal and spatial beam profile during propagation to 
          check for self steepening etc (read papers first!)
	
- improve nonlinear effects
	+ read through Luna documentation again ("Kerr_env_thg"...)
	+ read papers 
	+ email Luna developers?

- carry out parameter scans and see if filamentation effects (saturation, 
  second spectral peak, beam deformation) set in 
	=> find conditions to maximise conversion efficiency (and 
	   minimise pulse duration)

	+ vary pressure
	+ vary gas 
	+ vary beam energy 
	+ vary CEO phase 
	+ turn ionisation on/off





