GOAL: improve THG simulation
-----------------------------

Step 1: 

	- reproduce Josina's existing simulation results (familiarisation with code)
	- tidy up code where necessary

		+ ensure origin of all variable values are clear
		+ fix input variables for "Nonlinear.PlasmaCumtrapz" in "responses"
		
	- produce suitable output visualisation for easy comparison with measured data


Step 2:

	- set up "Fields.DataField" to represent input pulse (replacing "Fields.GaussGaussField")
		
		+ extract spectral phase from FROG data  (read through Jupyter Notebook)
		+ combine spectral phase with frequency and spectral power density
                  into a data file 