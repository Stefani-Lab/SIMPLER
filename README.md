# SIMPLER
 run_SIMPLER is a graphical user interface that runs in Matlab (version 2015a or later with Curve Fitting Tool installed). This app allows users of SIMPLER to perform all necessary operations to decode the axial positions of single molecules directly from 2D-SMLM-TIRF data. 

 The software also includes modules to perform the following operations: 

 - Determination of N_0 from 2D-SMLM-TIRF data of emitters bound or adsorbed to the coverslip
 - Calculation of the excitation intensity profile from 2D-SMLM-TIRF data of molecules spread all over the field of view; the list must include the background or offset information for each emitter.
 - Adjustment of the calibration parameters θ_i, α, and N_0 using the SIMPLER 3D reconstructions of standard structures of well-defined, known geometry as feedback.
 This software package includes several Matlab scripts and auxiliary functions, which implement the computational algorithms for the method SIMPLER described in the following publication:

 "Three-dimensional total-internal reflection fluorescence nanoscopy with nanometric axial resolution by photometric localization of single molecules."

 In the software documentation ("Supplementary Software Documentation.pdf"), users will find a detailed description of how to load files, set the different parameters, calculate z-coordinates of single molecules, and run the different operations available. 
 
 To use this software in GUI form, set the folder containing the software files (SIMPLER Supplementary Software.zip) as the current directory in Matlab and run:
 run_SIMPLER
 in Matlab's command window.


 

 



