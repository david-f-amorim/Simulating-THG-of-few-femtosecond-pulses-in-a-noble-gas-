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

finished handover presentation; wrote "manuscript_spectra.py" ; added new pressure scan plots (2d) ; 
adjusted code for relevant plot requirements; started on final presentation and report in earnest


TO DO
=====

- clean up files and directories [if significant changes: updated handover PP]
- make necessary changes to code [if significant changes: updated handover PP]
	+ write out UV timing data as file ?! 
	+ write code for Gaussian fit to UV timing data (TL duration)
        + arrange some plots with 2x2 subplots / increase fontsize [(chirp: Ar/Ne, spec/I(t)), (profile/I(t), ion/no-ion)]
	+ write code to compare densities (new, old, grad) 
        + change heatmap/contourf of spatiotemporal plots 
- process old chip COMSOL data 
- run again for plots 
	+ shift lam0 to 800nm 
        + produce all relevant plots (first for own use):
		* energy comparison (pres. scan): Neon and Argon 
                * spectra comparison: Neon and Argon
                * 2d comparison: Neon and Argon 
                * UV profile evolution: Argon 
                * self-steepening and temporal profile with and without ionisation: Argon 
                * chirp scan spectra: Argon and Neon 
                * gas density comparison: outside THG region and THG region-only
		* old-v-new comparison (spectra): Argon and Neon

- update presentation with plots and re-write where necessary 
- re-write report and update plots	

- then: do remaining plots for manuscript
	
