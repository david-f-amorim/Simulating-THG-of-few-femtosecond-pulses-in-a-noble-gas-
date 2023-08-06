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
2	Ar	100	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT		
3	Ar	150	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT 
4	Ar	125	 yes	f	0.0	0.1-5.1	   yes	yes	done	 SAT
5	Ar	175	 yes	f	0.0	0.1-4.4	   yes	yes	done	 SAT; DIES AT 4.5!
6	Ar	200	 yes	f	0.0	0.1-3.1	   yes	yes	done     SAT; DIES AT 3.2	 
7	Ar	225	 yes	f	0.0	0.1-2.7	   yes	yes	done     SAT; DIES AT 2.8	 
8	Ar	250	 yes	f	0.0	0.1-1.5	   yes	yes	done	 SAT; DIES AT 1.6 
9	Ar	275	 yes	f	0.0	0.1-1.7	   yes	yes	done	 SAT; DIES AT 1.8
10	Ar	300	 yes	f	0.0	0.1-1.4	   yes	yes	done	 SIGNS OF SAT; DIES AT 1.5 


4	Ar	75	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT 		
5	Ar	150	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT 
6	Ar	300	 no	f	0.0	0.1-5.1	   yes	yes	done	 SAT

7	Ne	75	 yes	f	0.0	0.1-35.1   yes	yes	done	 SAT 		
8	Ne	150	 yes	f	0.0	0.1-30.1   yes	yes	done	 SAT 
9	Ne	300	 yes	f	0.0	0.1-20.1   yes	yes	done	 SAT 


7	Ne	75	 yes	f	0.0	10.1-35.1   no	no	done	 SAT ; ADD TO EXISTING DIR; UPDATE EN&EF 		
8	Ne	150	 yes	f	0.0	10.1-30.1   no	no	done	 SAT ; ADD TO EXISTING DIR; UPDATE EN&EF 
10	Ne	200	 yes	f	0.0	0.1-20.1    no	no	maxw     SAT	 		
11	Ne	400	 yes	f	0.0	0.1-20.1    no	no	maxw	 SAT ; DIES AT 16.1


12	Ne	75	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT 
13	Ne	150	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT
14	Ne	300	 no	f	0.0	0.1-10.0   yes	yes	done	 NO SAT


		next steps: - try to get Neon to saturate!! [find peak with 0.1 accuracy]
			    - investigate CEO: 150mW Ar f ion, 0.1 to 5.1 step 0.5 ; phi vals: pi/4, pi/2, 3pi/4, pi, 5pi/4, 3pi/2, 7pi/4 ,2pi
			   

- plan analysis of pressure scan outputs
	+ for the different cases, look at beam at saturation point in more detail (and at breakdown point! also compare to non-ionised
	  case to see difference -> which effects are due to ionisation?)
		(e.g. beam shape; self-steepening: leading or trailing edge...)	
			=> re-write self-steepening code: show IR only! [asked Josina for wavelength range]
	+ see if relationship can be found between beam energy/CEO phase and saturation pressure/energy 
			- fit curves: beam energy vs peak pressure; beam energy vs peak energy; same for CEO phase
	+ study effect of ionisation (compare on/off) -> probably only use once interesting effect has been found in ionised case,
          to check if effect due to ionisation or not; also make large-scale comparison with existing non-ion data

	+ compare Neon to Argon (different behaviour during saturation?) 

	+ Thoughts so far:
		@ filamentation effects to appear, just at higher pressures and energies than expected 
		@ the higher pressure (~x4) can be explained by inaccurate pressure measurements 
		@ the higher energy is not unphysical but difficult to explain?


- improve density distribution model
	+ start looking at COMSOL [open file]

