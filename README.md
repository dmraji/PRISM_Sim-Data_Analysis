# PRISM_Sim-Data_Analysis
Matlab code for data analysis of PRISM_Sim Geant4 output files

This code is responsbile for processing outputs from Geant4 (in text, from tuples).

As of now, this program does not work with Octave.

The testing for this program has been primarily undertaken with Cs-137 in mind, 
  and has been primarily designed to process the output of 662keV particle gun
  Geant4 outputs as such. However, this is not to say that this program will only
  work with 662keV outputs. A few paramters may need to change, a few lines commented
  out. 
  Theoretically, this would not be necessary - electronic transport simulation,
  if undertaken, would provide for better understanding of trapping and other
  non-uniform effects. As of now, they are compensated by isotope-specific algorithms.
  Immediate future work should/will involve segreagating these.
  
Use with outputs from my Single-Det_Simulations Geant4 program or Dan Hellfeld's
  Imaging_Simulations (PRISM); Note that with the latter, small adjustments might
  be needed to read in data vectors - I am not sure at this time whether he has 
  approved my pull request to add in phi & theta angle readout and messenging.
  
Comments in this program are up-to-date.

