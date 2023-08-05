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
to run pressure scans; started data collection; started to implement 
scan visualisations;  


- read papers before meeting!


- visualise pressure scan results:
	+ add options to compare energy/efficiency output from different scans 
	+ fix spectrum plot to only show ~10-15 pressure plots at once 
	+ think about comparing spectra from different scans 

		
- carry out parameter scans to see if/when filamentation effects (saturation, 
  second spectral peak, beam deformation) set in 
	
	=> first generate lots of scan data for general analysis, 
	   then look at individual cases in more detail

	Gas	Energy   ion	kerr	phi	p_range	  Git?	Sync	Status	Notes
	------------------------------------------------------------------------------
1	Ar	75	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT 		
2	Ar	150	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT 
3	Ar	300	 yes	f	0.0	0.1-1.4	   yes	yes	done	 DIES AT 1.5 AND UPWARDS; SIGNS OF SAT
4	Ar	75	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT 		
5	Ar	150	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT 
6	Ar	300	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT

7	Ne	75	 yes	f	0.0	0.1-10.0   yes	yes	done	 NO SAT 		
8	Ne	150	 yes	f	0.0	0.1-10.0   yes	yes	done	 NO SAT 
9	Ne	300	 yes	f	0.0	0.1-10.0   yes	yes	done	 NO SAT 
10	Ne	75	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT 
11	Ne	150	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT
12	Ne	300	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT

13	Ar	100	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT		
14	Ar	125	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT
15	Ar	175	 yes	f	0.0	0.1-4.4	   yes	yes	maxw	 SAT; DIES AT 4.5!
16	Ar	200	 yes	f	0.0	0.1-5.1	   no	no	maxw	 INTERRUPTED
18	Ar	225	 yes	f	0.0	0.1-5.1	   no	no	maxw	 INTERRUPTED
19	Ar	250	 yes	f	0.0	0.1-5.1	   no	no	maxw	 INTERRUPTED
20	Ar	275	 yes	f	0.0	0.1-5.1	   no	no	maxw	 INTERRUPTED
		
	
		next steps: - go higher in pressure for Ne? [10.0-15.0 in steps of 1.0?]
			    - investigate more energies (maybe in steps of 25? or 50? go higher for Ne! also larger pressure steps!); 
			    - investigate CEO
			    => adopt much higher step size to speed up process? (0.1 to 5.1 in steps of 0.5?)

			LOOK INTO OTHER SERVER!

- plan analysis of pressure scan outputs
	+ for the different cases, look at beam at saturation point in more detail 
		(e.g. beam shape; self-steepening: leading or trailing edge...)	
			=> re-write self-steepening code: show IR only! [asked Josina for wavelength range]
	+ see if relationship can be found between beam energy/CEO phase and saturation pressure/energy 
	+ study effect of ionisation (compare on/off)

	+ Thoughts so far:
		@ filamentation effects to appear, just at higher pressures and energies than expected 
		@ the higher pressure (~x4) can be explained by inaccurate pressure measurements 
		@ the higher energy is not unphysical but difficult to explain?


- improve density distribution model
	+ start looking at COMSOL [open file]

