This is a brief guide on Fcyl_CCM.py

Fcyl_CCM.py was modified from Antdiatom_model.py to simulate CO2/HCO3- fluxes\
during photosynthesis, with updated parameters and constraints from experiments.\

It reads/calculate the parameters from the input file\
"Fc_3C_paras.csv" and a given column (default column nubmer= 1 (second column,\
 Python starts 0)), and simulates the carbon fluxes.
 #parameter_file = 'Fc_3C_paras.csv'
 #bpdict=get_parameter(parameter_file, col=1)

 One can specify a different column or change the parameters within the \
 Fc_3C_paras.csv file to see different simulation results under different conditions.
 When temperature is changed, one needs to make sure the kinetic and other parameters\
 reflect the temperature being used.

 The other files besides Fcyl_CCM.py and "Fc_3C_paras.csv" are the output files.
 Please note that the output files has concentration units in µM, \
 fluxes in mol/cell/s.

 The first version of the code was tested last on March 9th, 2022 within Python 3.7.6 using Spyder 4.0.1
 The updated version of the code was tested last on December 22nd, 2025 within Python 3.11.14 using Spyder 6.1.0

 Please Cite
 Li, Meng, and Jodi N. Young. "Extracellular carbonic anhydrase supports constitutive HCO3− Uptake in Fragilariopsis cylindrus regardless of temperature changes." Biorxiv (2022): 2022-09.

or
 Li, Meng, and Jodi N. Young. 2026, coming soon after peer review...

 if this code was helpful to your research.
