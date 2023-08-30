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
added option to norm plots; finished re-writing "scan_eval" 

Week 7
-------------------------------

finished handover presentation; wrote "manuscript_spectra.py" ; added now pressure scan plots (2d)



TO DO
=====

- find pressure scaling factors for both Argon and Neon: 2.5 ??
- adjust code for pressure scaling, where relevant [DONE??]
- cut p
- run simulations with old chip coms overnight 
- run Ne 400mW simulations for grad to higher pressures [use existing ones??!!]
- produce all relevant plots
	DONE: A,B,I for Ar and ne   
	      
	      missing: C,[D,E],F   [only show envelope, not carrie ->; done on maxwell ]
				 [add TL duration to chirp scan]
				

	running right now: Ar & Ne with and without ion pdf and png [covers D&E] ; Ar & Ne chirp scan [covers F]
	to do: C (wait COMSOL); 

	      G: scrapped
              H: on hold 


	RRE-RUN SIMULATIONS WITH 800nm ->> or just shift everything by 20nm ?? 

	REMOVE PULSE DURATION FROM PLOTS / PRINT OUT UV TIME DOMAIN FILE!!! ; 

		PROBLEM: pulse duration very inaccurate!!!

	NOTE: FOCUS FOR NOW ARE PRESENTATION AND REPORT; MANUSCRIPT LATER; FINISH DOING PLOTS [MAINLY EXPORTING AND WAITING FOR COMSOL]
		THEN WRITE REPORT AND PRESENTATION; THEN GO BACK TO MANUSCRIPT!!!


	NOTE: GENERAL BLUESHIFT OF MEASURED VS SIMULATED SPECTRA ->> SIMPLY BECAUSE WRONG CENTRAL WAVELENGTH??


		UPDATE THE PLOTS USED IN REPORT AND PP!!!


DO PROMISING SPECTRA PLOT??!!!

- run "final simulations" with good parameters [and PPT!]
	+ calculate ratios of Ar/Ne energies [just numerically at sat points; DOES NE EVEN SAT IN EXP??]
	+ 
	
- work on "documentation"
	+ write report
	+ prepare final presentation [use CFEL template; 12 minutes; formal]
	+ FIND PICTURES OF NEW AND OLD CHIPS!!!

			
			