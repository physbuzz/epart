Code to solve for the exact time evolution of a large system of rigid spheres undergoing elastic collision. 

The goal is to get something working and make some quick molecular dynamics videos. I'm not going to worry about the exact solution / exact event driven molecular dynamics code. 
It would be fun to scale it up to millions and billions of particles on a GPU, but I'll save that for later!

For a given particle at position $x_i$, we only have a collision if $\Delta x_{ij}\cdot \Delta v_{ij}\lt 0$, $\Delta x^2\gt r^2$ (which I should always enforce, but maybe an initial condition is given where the particles overlap), and $(\Delta x\cdot \Delta v)^2-\Delta v^2(\Delta x^2-r^2)\gt 0$. 

Let's divide into a grid of M by M cells. We evolve by time $\Delta t$, and delta t should be small enough that we usually have 0 collisions, rarely have 1 collision, and almost never have 2 collisions (which I won't try to handle properly!)





