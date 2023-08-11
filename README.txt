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

	
Week 4 [IN PROGRESS]
-------------------------------- 
		
added spatiotemporal plots ; changed to ADK ionisation from PPT ionisation to 
avoid DomainError issues; added beam focus position as parameter



- start testing Josina's COMSOl data 
	DOES NOT SEEM TO ACCEPT z<0 ???

- reproduce original plots with overlayed measured UV spectra

- look into different gases (:He, :Kr, :Xe, :Air, :N2, :H2, :O2, :CH4, :SF6, :N2O, :D2)
	started on Maxwell :Ar, :Ne, :Kr 
		[300muJ, 1.0bar]

- prepare question for conversation with Chris

- for next round of pressure scans:
	+ state 1/e^2 intensity [W/m^2] in addition to beam power (started on this)
		- FIX TITLE/PARAMS (PADDING ETC)
                - WRONG ORDER OF MAGNITUDE?? PEAK INTENSITY ??? CUBIC ???
	+ give ratio of saturation pressures and beam energies for Ar and Ne
	+ potentially add scan analysis for different beam focus values 
	=> wait until next round discussion with Vincent & Josina



