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
UV energy; started exhaustive new round of parameter scans; 


- reproduce original plots with overlayed measured UV spectra [realistically at higher pressures!]			
- prepare questions for conversation with Chris [Wed]
	* documentation on nonlinear responses
	* modal decomposition for radially symmetric free-space [mathcal{T} is Hankel transform? Bessel modes?]
	* phase output? seems to work with real fields tho


	- carry out new round of pressure scans

		+ varied variables: [140 total scans, i.e. 7000 simulations]
			* beam power (75mW, 150mW,200mW, 300mW, 400mW) [5]
			* pressure (0.1bar, 0.2bar,...,5.1bar) [50]
			* gas (Ne, Ar, Xe, Kr, He, N2, N2O) [7]
                        * ionisation on/off [2]
			* COMSOL data/gradient data [2] (maybe only with Ar/Ne for now?)
		+ in scan_eval: 
			* replace phi scan with gas scan 
			* make distinction between coms and grad
	

		CARRY OUT LAST TESTS AND THEN START RUNNING SCANS!

		started: [20 rows]
		    1	- 075mW; ion; com;    Ar, Ne, Kr, Xe, He, N2, N2O    [all done; upload to Sync&Share + GitHub]
		    2	- 150mW; ion; com;    Ar, Ne, Kr, Xe, He, N2, N2O    [all done; upload to Sync&Share + GitHub]
   		    3	- 200mW; ion; com;    Ar, Ne, Kr, Xe, He, N2, N2O    [all but Xe  ]
		    4	- 300mW; ion; com;    Ar, Ne, Kr, Xe, He, N2, N2O    [all but Xe   ]
		    5	- 400mW; ion; com;    Ar, Ne, Kr, Xe, He, N2, N2O    [done: He, N2, Ne]
		
		    6   - 075mW; ion; grad;    Ar, Ne, Kr, Xe, He, N2, N2O   [done: Ar,He, N2, Ne ]  
		    7   - 150mW; ion; grad;    Ar, Ne, Kr,                              
                    8   - 200mW; ion; grad;    Ar, Ne, Kr, Xe, He, N2, N2O   
		    9   - 300mW; ion; grad; 
		    10  - 400mW; ion; grad;        

		  
		  [17 sessions running!]


- continue writing report 
	NOTE: startlight not actually part of beamline!!

	
	


