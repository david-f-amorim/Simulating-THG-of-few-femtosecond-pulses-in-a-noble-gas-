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
cleaned up "file_prepare.py"; finished gas scans; changed from on-axis to total intensity/energy; 
started re-writing "scan_eval.py"; now integrating along frequencies for intensity distributions; changed plot design;



- run "final simulations" with good parameters [and PPT!]
	+ calculate ratios of Ar/Ne energies [just numerically at sat points]
	+ work on plot list 

	
	
- work on "documentation"
	+ write report
	+ prepare final presentation [use CFEL template; 12 minutes; formal]
	+ prepare handover PP 
			[clean up scan_eval] ; upload gas_scans etc to Sync&Share (fixed param files)
			-> add functions for multi_gas comparison
			-> add quick notes section for control
			-> write PP slide

			go through THG_main again
			-> change some plot labelling
			-> if making changes for manuscript plots, make sure to update handover


