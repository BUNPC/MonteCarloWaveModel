# MonteCarloWaveModel
This code integrates Monte Carlo simulation with wave model to simulate time domain diffuse correlation spectroscopy. 
The Monte carlo module is similar to our previous Monte carlo code. The source code is: tMCimgLOT_ScatDir_TDDCS.c

The wave module is:
BaselineV7_AnalyzeData_multipleWavelength_twoLayers.m
This code requires the .his file from Monte Carlo simulations and an input IRF profile (we have included FWHM_NearGaussian_320ps.mat, which is the IRF profile we have measured experimentally).
