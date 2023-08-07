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

added different options for nonlinear model (use "ff" for now) ; 
added visualisation of pulse self-steepening; added option to overlay
input pulse spectrum; changed cell length to 3mm ; implemented option 
to run pressure scans; started data collection; implemented 
scan visualisations;  


- read papers before meeting!

- extend grid: from 100nm!

- interpret scan analysis!
	
- plan analysis of pressure scan outputs
	+ for the different cases, look at beam at saturation point in more detail (and at breakdown point! also compare to non-ionised
	  case to see difference -> which effects are due to ionisation?)
		(e.g. beam shape; self-steepening: leading or trailing edge...)	
		
	
	+ compare Neon to Argon (different behaviour during saturation?) 

	+ Thoughts so far:
		@ filamentation [don't use this term!] effects do appear, just at higher pressures and energies than expected 
		@ the higher pressure (~x4) can be explained by inaccurate pressure measurements 
		@ the higher energy is not unphysical but difficult to explain?


- improve density distribution model
	+ start looking at COMSOL [open file]

