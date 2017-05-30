# S Durbridge 
s.durbridge1@unimail.derby.ac.uk
Independent Engineering Scholarship 2017
7EJ998

Efficient Acoustic Modelling of Large Spaces Using Time Domain Methods

README

This CD contains two folders: Docs and Matlab.

The Docs folder contains the source code for the disserataion (written in LaTeX), a library of supporting literature, a bibtex library, graphics and so on.

The Matlab folder contains the reduced Matlab code for the study i.e. the final code of interest for this study. 
The Matlab folder contains five main subfolders, and some scripts. The subfolders contains the FDTD simulation material, the PSTD simulation material, the SFDTD simulation material, an Image Source modelling program and an mls generating package.

Each of the FDTD, PSTD and SFDTD folders contains a set of functions and scripts. The functions execute the main update equations, and some surrounding functionality for the execution of each program.
The [..]TD_Dtesting.mat files are those used for validating each method up to 5kHz. Both the 2D and 3D files are provided.
The [..]TD_Dtestingexec.mat files are those used for testing the execution speed of each method with different domain sized up to 500Hz.
These script should simply run in the matab environment using an apropriate command, or the run button. 

The scripts are devided into three main sections, initialization, execution and postprocessing. 
Initialization creates all the variables for the simulation.
Execution loops the algorithm over the number of time steps.
Postprocessing normalises the data and displays the time and frequency domain content.
'Runtime' or in-the-loop plotting for some scripts will be comented out for the sake of saving time while running the code. If this is of intrest, 

The MLS package conatined in the mls folder is references in the main literature as the work of Mark Thomas, and contains both MLS generation and analysis functions and a supporting document that intoduces the theory behind MLS.

The IMS folder contains an Image Source Method program written in Matlab, under the supervision of Dr Adam Hill, during the undertaking of the MSc in Audio Engineering. This folder also contains impulse response data for the domain used for validating the three numerical methods in this study.

The matlab scripts inside the matlab folder was used for organising and plotting some of the data presented in the report.
