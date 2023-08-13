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


RESPOND TO JOSINA'S EMAIL!


- start testing Josina's COMSOl data [new and old] 
	DOES NOT SEEM TO ACCEPT z<0, shift data!

- reproduce original plots with overlayed measured UV spectra [realistically at higher pressures!]

- look into different gases (:He, :Kr, :Xe, :N2, :H2, :O2, :CH4, :SF6, :N2O, :D2) [300muJ, 1.0bar]
		
		to compare: -plot all UV spectra on top of each other 
			    - plot "bar charts" with pulse duratons and 
				UV energies
				=> also plot IR beam depletion as measure of ion.?
			put other plots on PP?	

- prepare questions for conversation with Chris [Wed]



	
	


