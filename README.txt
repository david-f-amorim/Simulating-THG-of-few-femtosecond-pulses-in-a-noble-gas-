DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; carry out parameter scans to find optimal 
	  input beam conditions for desired UV output beam;   
		 

Week 1
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

Week 2
----------------------------------

provided option to base input pulse on measured data ; provided 
option to feed in custom gas density profile (untested!) ; 
provided option to overlay measured UV spectrum 
	

Week 3 & 4
--------------------------------  

added different options for nonlinear model (use "ff" for now) ; 
added visualisation of pulse self-steepening; added option to overlay
input pulse spectrum; changed cell length to 3mm ; implemented option 
to run pressure scans; started data collection; implemented 
scan visualisations;  extended wavelength grid to start at 100nm (from 200nm); 
analysed scan results; added spatiotemporal plots ;  



- reproduce Reiter plots [see old code]
	-> RUN WITH/WITHOUT IONISATION AND COMPARE
	-> PUT ON PP!



- run scans with different focus values
	ADDED PARAMETER [ANNOTATE ON PLOTS?]

MOSTLY LABELLING ISSUE?
- give saturation pressures and beam energies in ratios of Ar/Ne 
- change beam power to intensity 


DONE!

- list points for discussion w Chris (Luna)
		+ documentation for nonlinear functions
		+ phase information
- look into spectral and temporal input beam profile disagreeing [compare Luna, Speck.dat, experimental] 
		


	
Week 4 [IN PROGRESS]
-------------------------------- 		
	


- improve density distribution model
	+ start looking at COMSOL [open file]

