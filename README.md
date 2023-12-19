# Diffusiophoresis
MATLAB code for simulating the diffusiophoretic effect of charged ions in a fluid with a concentration gradient

diffusion.m - Nernst-Planck equation solver, transient simulation for calculating ion concentration values. Simultaneously solves for potential at each point in time as well. 

initialPotential.m - plots initial potential field for the initial conditions given in code "diffusion.m". (essentially just the first time step in "diffusion.m")

fullyDeveloped.m - yields steady-state and fully-developed solution for ion concentrations. 

newFluxGridDiffusion.m - Nernst-Planck equation solver with a different grid for the mesh.

new_approach.m - new solver for calculating ion concentration values using a different method than "diffusion.m".

new_approach_myCase_fullyDeveloped.m - steady-state and fully-developed solution for ion concentrations using the different method.
