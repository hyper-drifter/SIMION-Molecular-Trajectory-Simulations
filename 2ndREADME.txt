Taylor Grubbs 3/12/16

MaxPlanckQuad_v2 uses a certain stark slope defined by the state of the ND3 molecule.
In the MPI paper, these states were defined by 4 different shifts were simulated: 1.2, .6, -1.2, and -.6 wavenumbers at 100 kV/cm.

In SIMION units this converts to 1.44*10^-8 eV/(V/mm) for 1.2 wavenumbers and 7.2*10^-9 eV/(V/mm) for .6  wavenumbers

Energy calculations were made assuming a linear Stark shift defined by W = sE where s is the slope of the stark shift and E is the electric field.

This version runs the same particles through the simulation at different frequencies and records the resulting data from each run.
