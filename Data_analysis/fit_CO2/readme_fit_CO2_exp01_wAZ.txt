This is a brief guide on fit_CO2_exp01_wAZ.py for analyzing CCM experimental data\
acquired with AZ, eCA inhibitor.

The analysis depend on the input file, sample data was given:
"Fcyl_3C_S40_exp01_c.csv"

This guide covers the first of the two experiments: \
1) CCM experiment with AZ; \
2) CCM experiment without AZ.
If eCA is not a concern, the first experiment is sufficient to ask\
questions about CO2/HCO3- usage by a diatom or other phytoplanktons.

By default, the scripts analyze the sample data. It is recommended to analyze the\
data section by section as the experiment goes:

The general steps are:
 1, identify an experimental section,
 2, do the analysis/calculations of the section
 3, store the values useful for next section
 4, go to step 1 and repeat the process until the end of the experiment.

================== fit_CO2_exp01_wAZ.py Section 1 ==============================
adding 13C, 18O - labeled NaHCO3 to a given buffer, the dehydration of HCO3-\
and 18-O exchange between H2O and HCO3- can be used to simulate the 13C-CO2 \
(18O and 16O) signal and curve fit the signal to certain k1, k2 values so that
    CO2--k1--> HCO3-     and CO2  <---k2---HCO3-

#### use 49CO2 to determine time of addition of NaHCO3
One needs to make sure time variables, t0 and t1, cut the section of \
passive HCO3- dehydration and 18O exchange (before adding cells).
Then,
#### HCO3- dehydration simulation with approximate k1, k2 values and fit
The following lines calculate passive process parameters:
    k1,K_bulk, DIC0, f180 = fit_passive_3(passive.iloc[:])
    k2 = k1/K_bulk

Alternatively, one can use the fit_passive function to calculate k1, k2...


================== fit_CO2_exp01_wAZ.py Section 2 ==============================
Adding cells to reaction chamber leading to faster 18O-exchange
#### Add cells, fc, fb, kcf... ####
There are a few parameters that need to be experimentally specific
    N = 6.04e8*0.020    #total number of cells
    Ve = 1226*1e-3      #rxn cuvette volume in mL
    Vc = 41.6e-12        #cellular volume in mL, measured average volume at 3°C
Also, cell diameter is used to calculate Vs, surface volume.
#measured average effective diameter ~4.30µm @ 3°C

This section estimates time points of events, that is adding cells, light on,\
light off, light on, light off...
those timepoints are stored in index list variable t.
One needs to make sure those time points cut the sections of interest, for instance,
dark0 is the data section from adding cells to first light on, which is used to\
estimate:
    kcf, Kcyto, fc, fb = fit_dark_1C(dark0, y01)
    kcr = kcf/Kcyto
with fc, fb and cell size (diameter), fc_M and fb_M can be calculated.


================== fit_CO2_exp01_wAZ.py Section 3 & 4 ==========================
Section 3
#### first light period
Assigning light1, the first light on period.
Factoring in MIMS consumption residual after Ar correction.
Calculate HCO3- before turning light on:
    HCO3_dark0 for first light period, calculated using calc_DIC_dark function.
Then calculate PS, DIC usage using calc_PS_Urates function.

Section 4
#### calculate light2 rates
assign light2, the second light on period
Calculate HCO3- before turning light on:
    t4_HCO3 for second light period, estimated assuming ∆O2 = ∆CO2+∆HCO3-.
Then calculate PS, DIC usage using calc_PS_Urates function again.

Steady state values were averaged for the ending 5 minute during a light period

=======Output=======
files:
1, "*wAZ.pdf":
Plots of PS rate, DIC usage calculation and results.
2, "*1st_light_wAZ.csv" and 3, "*2nd_light_wAZ.csv":
data storing smoothed CO2/HCO3- O2 data and Uc, Ub...

The fit_CO2_exp01_wAZ.py was tested on March 10th, 2022, within Phython 3.7.6 using\
Spyder 4.0.1

Known issue:
SettingWithCopyWarning:
A value is trying to be set on a copy of a slice from a DataFrame.
Try using .loc[row_indexer,col_indexer] = value instead
See the caveats in the documentation: https://pandas.pydata.org/pandas-docs/stable/user_guide/indexing.html#returning-a-view-versus-a-copy
