# xlzd_reqtask_nufloor
Repository for the Requirements Taskforce neutrino fog work (2025)

Group Members: 
 * Dan Tovey
 * Robert James
 * Maike Doerenkamp
 * Knut Morå

## Contents: 

  * template_generation # Folder containing scripts to generate expected observable distributions
  * data #to store spectra, templates (<10-ish Mb)
  * likelihoods #store likelihood configurations of all runs
  * run_configurations #store run configurations of the statistics runs
  * results #store (abbreviated) run results here


## Run versions
Each likelihood, run config and results should be associated with a version name to allow easier organisation, listed below
  * k0.0: A very simplified ER+nu+WIMP likelihood using templates from Rob
  * r0.2: Rob ran discovery significance projections and defined a "good" and "bad" detector state. General detector config found in https://github.com/FlamTeam/flamedisx/blob/RJ-XLZD_simple/flamedisx/xlzd/xlzd.py :  
    * Good: 80V/cm drift, 7.5 keV gas field, 10ms electron lifetime, 0.27 PMT quantum eff. 
    * Bad: 25 V/cm drift field, 6 keV/cm gas field. 
  * 0.5: first unified run+template generation to be set and varied for the XLZD collaboration meeting
  * 0.51: parameters for example XLZD plot

## Use
You can install the fogtask package as another python package-- clone this repo, and run ```pip install -e .```in this folder. 

To generate templates, you can then run: 


```
from fogtask.fogtask import *
generate_all_wimp_templates(version='v0.5',
                            n_samples = int(1e7), 
                            file_name_pattern="templates_{version:s}_{parameter_string:s}")
```

