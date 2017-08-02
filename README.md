
 ##   WindFarmSimulator (WFSim)
Developed by Boersma et al., Delft University of Technology, 2017
            


## Summary:
WindFarmSimulator (WFSim) is a medium-fidelity, control-oriented wind farm model based on the two-dimensional Navier-Stokes equations. It is currently actively developed at the Delft University of Technology by Sjoerd Boersma.

## Quick use:
Open WFSim.m with any recent version of MATLAB. Follow the instructions therein to perform simple simulations of various wind farm scenarios. Note: if you are missing files such as 'system_input.mat', make sure you have downloaded the relevant files in the data_SOWFA folder. Specifically, for the default simulation scenario, download the files in the YawCase3 folder.
	
## Folder hierarchy:

	/bin/:          contains all the functions and scripts used by WFSim.
	/bin/analysis/: files not essential to WFSim, but used for debugging and validation.
	/bin/archive/:  outdated files no longer used by any script.
	/bin/core/:     files essential to the workings of WFSim, required for any simulation.
	
	/data_PALM/:    high-fidelity simulation data from Hannover's LES code 'PALM'. This can be used for model validation.
	/data_SOWFA/:   high-fidelity simulation data from NREL's LES code 'SOWFA'.    This can be used for model validation.
	/data_WFSim/:   simulation data from WindFarmSimulator itself, useful for development, debugging and verification.
    
	/documentation/: literature on the technical details of WFSim, such as the derivation and the turbulence model.
		          Please keep in mind that more recent literature may be available online.
	/libraries/:    external libraries used in WFSim. All copyright goes to the respective authors.
	
## Debugging:
For any serious issues, reach out to us on the Github page. 

All credit goes to the Delft University of Technology. WFSim was written by ir. Sjoerd Boersma with support from ir. Bart Doekemeijer and under the supervision of dr.ir. Jan-Willem van Wingerden.             
