DESY Summer Programme 2023
============================

GOAL:     simulate THG of a few-femtosecond IR laser beam in a 
----      noble gas; 
	  carry out parameter scans to find optimal input beam 
          conditions for desired UV output beam;   
		 


Stage 1: DONE [Week 1-2 ]
--------------------------
	
did background reading ; cleaned up existing code ; 
provided visualisation functions ; provided output handling 

	(MAYBE: reproduce Josina's results ?) 


Stage 2: IN PROGRESS [Week 2 - ?]
----------------------------------

	- set up "Fields.DataField" to represent input pulse (replacing "Fields.GaussGaussField")
		
		+ current problem: GaussGaussField is type "SpatioTemporalField", DataField is type "DataField" ("TimeField")
			-> function "doinputs_fs!" (executed as part of "setup") [line 235 of Luna.jl] requires the input field to be a 
			   "SpatioTemporalField"

		+ APPROACH 1: USE A DIFFERENT VERSION OF "setup"??

			- there are 8 different versions of "setup" in total; of which 4 work with grid::Grid.RealGrid [as opposed to 			          grid::Grid.EnvGrid] (changing grids does not make a difference!)
			
				A. [Luna.jl line 141]: uses "inputs" in "doinput_sm!"

					function setup(grid::Grid.RealGrid, densityfun, responses, inputs, βfun!, aeff;
              						 norm! = NonlinearRHS.norm_mode_average(grid, βfun!, aeff))

					- aeff : function that gives effective area of mode at position z ?? [see Modes.Aeff]
							Modes.Aeff : requires m::AbstractMode and z 
					- βfun!: function that reads in a zero-array fills it with the β coefficient at z?? [see 							 Modes.β]
							 Modes.β : requires m::AbstractMode, omega, and z 
								->> which omega???
								->> based on Modes.neff which is EMPTY??? / maybe reference to 
								    Capillary ?

					- both of these only relevant for feeding into "NonlinearRHS.norm_mode_average" !

					- HOW TO SET UP Abstract:Mode --> ALL THESE MODES ARE FIBRE MODES! WILL NOT WORK FOR FREESPACE!!

	 						I.  StepIndexMode [StepIndexFibre.jl] (seems unusable for freespace)
							II. SimpleMode [SimpleFibre.jl] (requires Taylor expansion in beta..)
								-> probably not usable for freespace ??
							III. RectMode [RectModes.jl] (based on fibre props. -> not for freespace)
							IV.  ZeisbergerMode [Antiresonant.jl] (not for freespace ??)
							V.   MarcatiliMode [Capillary.jl] (not for freespace ??)
				
				B. [Luna.jl line 191]: uses "inputs" in "doinput_mm!"

					function setup(grid::Grid.RealGrid, densityfun, responses, inputs,
               						modes::Modes.ModeCollection, components;
              						full=false, norm! = NonlinearRHS.norm_modal(grid))

					->> SAME PROBLEM AS WITH A: modes cannot be specified !!!

				C. [Luna.jl line 247]: uses "inputs" in "doinputs_fs!" (currently used!)

					function setup(grid::Grid.RealGrid, q::Hankel.QDHT,
               						densityfun, normfun, responses, inputs)

			{	D. [Luna.jl line 281]: uses "inputs" in "doinputs_fs!"
	
					function setup(grid::Grid.RealGrid, xygrid::Grid.FreeGrid,	
							densityfun, normfun, responses, inputs)               } 
			
			- version C is currently used, problem with "doinputs_fs!"; same problem with version D so no 
			  point in changing to D 

				+ "doinputs_fs!" requires "inputs::Fields.SpatioTemporalField"

				-> sticking to current Code requires changing input to SpatioTemporalField

			- versions A/B (i.e. "doinputs_sm" and "doinputs_mm") allow for "inputs::Fields.TimeField"



		+ APPROACH 2: BRUTE FORCE DATA INTO "SpatioTemporalField" !!!

			->> requires carrier-envelope offset phase (CEO) [scalar, not array] ???

				set to zero for now ; put as input parameter maybe 
			
				use os.walk to search through all data files for "CEO" or "CEP" ??
			
			->> requires temporal shift from grid time zero (=0??)
			->> requires propagation distance from the waist,  "propz"			


			->> work with temporal intensity  from FROG data [Ek ?!]

				-> what about the spatial dependence ??
				-> assume Gaussian in space ??!! [multiply temporal part by Gaussian in space]

			->> GaussGaussField seems to be only instance 
			
			->> IDEA: use c spline interpolation [Maths.CSpline or Bspline ??] to generate callable Ishape (?)
				
					->> use envelope only -> hilbert transform of E(t)
		

	- compare input and output with measured data 

		+ NIR spectrum at input compared to measured data
		+ UV spectrum at output compared to measured data
		+ ....


Stage 3: TO DO
-----------------------------

	- vary pressure, beam intensity, etc. to study effect on THG efficiency [Luna parameter scans...]
		-> try to reproduce Josina's measured pressure scans (saturation, pulse shape change, second spectral peak)

	- set up own density function ? set up non-Gaussian spatial profile ?