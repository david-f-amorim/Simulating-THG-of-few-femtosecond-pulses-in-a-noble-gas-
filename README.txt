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


- start testing Josina's COMSOl data 
- start investigating different beam focus values (propz)


	
Week 4 [IN PROGRESS]
-------------------------------- 		
	
- for next round of pressure scans:
	+ state 1/e^2 intensity [W/m^2] in addition to beam power (started on this)
	+ give ratio of saturation pressures and beam energies for Ar and Ne
	+ potentially add scan analysis for different beam focus values 
	=> wait until next round discussion with Vincent & Josina

- improve density distribution model
	+ start looking at COMSOL [open file]

