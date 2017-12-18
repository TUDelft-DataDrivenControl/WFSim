
 ##   WindFarmSimulator (WFSim)
Developed by Boersma et al., Delft University of Technology, 2017
            


## Summary:
WindFarmSimulator (WFSim) is a medium-fidelity, control-oriented wind farm model based on the two-dimensional Navier-Stokes equations. It is currently actively developed at the Delft University of Technology by Sjoerd Boersma and Bart Doekemeijer. The most recent publication on WFSim can be found here: https://www.wind-energ-sci-discuss.net/wes-2017-44/

## Quick use:
Open WFSim.m with any recent version of MATLAB. Follow the instructions therein to perform simple simulations of various wind farm scenarios. Missing files will be downloaded automatically on first run, so make sure you are connected to the internet.
	
## Folder hierarchy:

	/bin/:          contains all the functions and scripts used by WFSim.
	/bin/analysis/: files not essential to WFSim, but used for debugging and validation.
	/bin/archive/:  outdated files no longer used by any script.
	/bin/core/:     files essential to the workings of WFSim, required for any simulation.
	/data_LES/:     high-fidelity simulation data from Hannover's LES code 'PALM' and NREL's LES code 'SOWFA'. 
                    These datasets will be downloaded auotomatically by meshing.m and can be used for model validation.

	/documentation/: literature on the technical details of WFSim, such as the derivation and the turbulence model.
		             More recent literature is available online: https://www.wind-energ-sci-discuss.net/wes-2017-44/

	/libraries/:    external libraries used in WFSim. All copyright goes to the respective authors.
	
## Debugging:
For any serious issues, reach out to us on the Github page. 

All credit goes to the Delft University of Technology. WFSim was written by ir. Sjoerd Boersma with support from ir. Bart Doekemeijer and under the supervision of dr.ir. Jan-Willem van Wingerden.             
