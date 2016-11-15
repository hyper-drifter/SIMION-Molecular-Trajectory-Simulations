# SIMION-Molecular-Trajectory-Simulations
In the Spring semester at UNC Greensboro I worked on a project involving the trapping of dipolar molecules via time-varying electric fields (the Stark effect). The simulation software that I used was the SIMION ion optics program. Here I have posted the lua files used to build and run the simulations as well as the Mathematica notebooks used to analyze the resulting data. The first goal of the project was to recreate the results found in the following paper https://arxiv.org/pdf/physics/0310046v2.pdf. The results that I obtained were remarkably similar to those reported in the above paper.

Usage: 

(1) After starting SIMION, run the files titled "maxPlanckFirstStageBuild_v2.lua" and "maxPlanckSecondStageBuild_v2.lua." These will build the necessary PAs needed to simulate the Stark effect on the ND3 molecule. 

(2) Open the "MaxPlanckQuad_v2AR.iob" file to open the workbench.

(3) Fly'm! Try viewing from other perspectives to see the micro/macro-motion of the trajectories.
