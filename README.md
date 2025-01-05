# Bachelor Thesis Philip Mundt

This repository contains the scripts that allow to execute the ftINIT and sMOMENT methods as a workflow in order to create tissue specific models from the Human1 model. 

In order to use the ftINIT & sMOMENT methods, please clone the Human-GEM & autopacmen repositories parallel to this repository. 
The scripts in the subfolder "sMOMENT/updated_autopacmen_scripts" have to be copied into autopacmen/autopacmen/submodules to replace the scripts there.

The packages needed for the execution of all python scripts are contained in the context_GEM conda environment. It can be created with the "context_GEM.yml" using 
conda env create -f context_GEM.yml 

As a starting point, the "workflow_setup.ipynb" file should be executed to create folders for input data.
