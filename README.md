# SIMPLER
 run_SIMPLER is a graphical user interface that runs in Matlab (version 2015a or later with Curve Fitting Tool installed). This app allows users of SIMPLER to perform all necessary operations to decode the axial positions of single molecules directly from 2D-SMLM-TIRF data. There is also the possibility to download SIMPLER as a standalone application, for users who do not have a Matlab license (URL: ).

 The software also includes modules to perform the following operations: 

 - Determination of N_0 from 2D-SMLM-TIRF data of emitters bound or adsorbed to the coverslip
 - Calculation of the excitation intensity profile from 2D-SMLM-TIRF data of molecules spread all over the field of view; the list must include the background or offset information for each emitter.
 - Adjustment of the calibration parameters θ_i, α, and N_0 using the SIMPLER 3D reconstructions of standard structures of well-defined, known geometry as feedback.

 This software package includes several Matlab scripts and auxiliary functions, which implement the computational algorithms for the method SIMPLER described in the following publication:

 "Three-dimensional total-internal reflection fluorescence nanoscopy with nanometric axial resolution by photometric localization of single molecules."

 In the software documentation ("Supplementary Software Documentation.pdf"), users will find a detailed description of how to load files, set the different parameters, calculate z-coordinates of single molecules, and run the different operations available. 

# Requirements
 Matlab 2015a and newer, Curve Fitting Tool installed. 
 
 For users who do not have Matlab license, Matlab Runtime 9.6 (R2019a) is requiered, which is free and can be downloaded from the web (https://se.mathworks.com/products/compiler/matlab-runtime.html). Alternatively, if users directly download the standalone SIMPLER application, the installer will detect whether Matlab Runtime is installed or not, and download it if it cannot find it.

# Usage
 
 For users without Matlab license: 
 
 The standalone application installer will create a folder called "application", which will contain the "run_SIMPLER.exe" executable. By clicking on this file, SIMPLER will automatically run.
 

 For users with Matlab license:


 To use this software in GUI form through Matlab, set the folder containing the software files (SIMPLER Supplementary Software.zip) as the current directory in Matlab and run in Matlab's command window:

 run_SIMPLER

The "Example data" folder includes raw data of microtubules, spectrin rings and nuclear pore complex, imaged in two different microscopes. The parameters needed for the 3D reconstruction of each example can be found in "Parameters.xlsx".

# How to cite
 If you used the code/program for your paper, please cite

# Contact
 Feel free to contact us with comments or suggestions. Use any part of the code that suits your needs.

 alanszalai@gmail.com
 
 a.szalai@cibion.conicet.gov.ar
