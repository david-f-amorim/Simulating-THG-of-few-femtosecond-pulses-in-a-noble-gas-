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
	

Week 3
--------------------------------  

added different options for nonlinear model (use "ff" for now) ; 
added visualisation of pulse self-steepening; added option to overlay
input pulse spectrum; changed cell length to 3mm ; implemented option 
to run pressure scans; started data collection; implemented 
scan visualisations;  extended wavelength grid to start at 100nm (from 200nm); 
analysed scan results;  

	
Week 4 & 5
-------------------------------- 
		
added spatiotemporal plots ; changed to ADK ionisation from PPT ionisation to 
avoid DomainError issues; added beam focus position as parameter; changed 
fundamental wavelength to 730nm; added peak intensity information to pressure scans ; 
prepared table of relative Ne/Ar saturation pressures and energies ; 
compared performance of different gases; extended scan_eval to also investigate pulse duration; 
added option to feed in COMSOL data; extended scan_eval to also investigate position of peak
UV energy; started exhaustive new round of parameter scans; compared COMSOL density profile to grad profile; compared density profiles of different gases; extended scan_eval to compare different gases as well as COMSOL vs grad; 


Week 6
--------------------------------
extendend previous pressure scans to O2 and SF6; implemented option to add chirp to the input pulse; 
carried out minor chirp scan; added option to toggle between PPT and ADK; started working on handover PP; 
cleaned up "file_prepare.py"; finished gas scans


- run "final simulations" with good parameters [and PPT!], based on measurements 
	+ calculate ratios of Ar/Ne energies 
	+ reproduce measured data 
	+ show what is "happening inside"; maybe turn ionisation on/off and compare... ; self-steepening...
	+ play with chirp/propz/CEP/pressure/beam power/beam waist...?
	+ WAIT FOR NEW COMSOL DATA!!

- work on "documentation"
	+ write report
	+ prepare final presentation 
	+ prepare handover PP [clean up scan_eval]
			-> change function and variable names!!
			-> rewrite!!! make more general!!





	
	


