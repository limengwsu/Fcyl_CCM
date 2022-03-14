# Fcyl_CCM
CCM modeling and MIMS data analysis for Fragilariopsis cylindrus
The Fcyl_CCM repository contains Python scripts for modeling _Fragilariopsis cylindrus_ CCM (carbon concentrating mechanism), as well as MIMS data analysis tools when following similar experimental procedures to investigate inorganic carbon usage by a diatom or other aquatic phototrophs.
## Software requirement
Python 3.7, Text editor Spyder 4.0.1 or above recommended. The Python scripts have been tested last time in Python 3.7.6 using Spyder 4.0.1. 

## Input and Output
This is folder dependent. Data_analysis is used for analyzing MIMS data from CCM experiments, while Kranz_2014 and Li_2022 have original _Fragilariopsis cylindrus_ CCM model and modified one. **Each of the following folders has its own readme_*.txt file**
### Data_analysis\

Please follow Gas_calibration then fit_CO2 instructions sequentially, each folder has its own __readme_*.txt file_
#### **Gas_calibration**\

Input: MIMS mass spec reading data in .asc format\

Output: csv files including calibration parameters and two experimental data in gas concentrations over time
#### **fit_CO2**\
Input: experimental data files generated from Gas_calibration\

Output: CCM calculation results and plots in csv and pdf format respectively.

### Kranz_2014
Input: parameter files stored as **Antdiatom_model_MATLAB/Antdiatom_large_S1.par**\

Output: Simulation results with plot in csv or pdf format

### Li_2022
Input: parameter file "Fc_3C_paras.csv"\

Output: Simulation results with plot in csv or pdf format

## The basic design
For data analysis, steps were coded into Python functions. To estimate/calculate parameters of interest, the "modeling plus curve fitting" approach was used.

For modeling, Python scripts follow original MATLAB code design published by Kranz et al 2014 (see reference). ![Dr. Brian Hopkinson](https://github.com/bmhopkinson) at University of Georgia, to his credit, created the original model and the predecessor ![Chloroplast Pump model](https://github.com/bmhopkinson/Ci_physiology_modeling/tree/master/Chloroplast_pump)
The original code was posted here with his permission.

## References
Please cite the following work if the CCM model is helpful:\

Kranz, Sven A., Jodi N. Young, Brian M. Hopkinson, Johanna AL Goldman, Philippe D. Tortell, and François MM Morel. "Low temperature reduces the energetic requirement for the CO 2 concentrating mechanism in diatoms." New Phytologist 205, no. 1 (2015): 192-201.\

Or cite the to-be-published work by Li M et al 2022 if the data analysis tools are useful for your work.\

Meng Li, Brian M. Hopkinson, Jodi N. Young 2022, under preparation...
