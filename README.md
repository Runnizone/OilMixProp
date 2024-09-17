# OilMixProp

OilMixProp 1.0 Manual by Xiaoxian Yang - July 18, 2024.
Feel free to contact Xiaoxian for any help.

1. For fluid property calculation
    Run FluidCalc.m
	See examples and comments in FluidCalc.m for more details

2. Add your own oil: TheOilName
   2.1 Prepare the experimental data in the folder \ExpData\TheOilName, using the template given in \ExpData\PAG68
       At least density data should be avaiable. 
   2.2 Add your oil case in Fitter_PureOil.m and run this code. Please see the comments inside.
	   The output is \ExpData\TheOilName\Fluid_Constants_TheOilName.txt
	   Another output is \Classes\Fluid_Constants_Fitted.txt
   2.3 Carry out fluid property calculation using FluidCalc.m

3. (not ready yet) Fit binary interaction parameter  (BIP) for an oil + refrigerant mixture
   Not very robust yet, if there are not experimental vaper pressure data, all BIPs are defaulted as zero.
   3.1 Prepare experimental data using the template given in \ExpData\Emkarate RL32_R1233zde
   3.2 Add the case in Fitter_CKij.m and run this code. Please see the comments inside.
   3.3 The output is \Classes\Bin_kij_fit.txt

4. There are also some functions to plot phase diagrams available. 
   4.1 PhaseDiag_pT.m (not quite robust...)
   4.2 PhaseDiag_px.m
   4.3 PhaseDiag_Tx.m
   They should work in most of the cases. Feel free to modify them for your best output. 

