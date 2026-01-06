Various n-dimensional particle simulation codes. (Note: Many features are only 
available in 2D.)

Originally this was just `epart`, intended to simulate nearly exact particle dynamics of particles that move with constant velocity between elastic collisions. `epart` is designed for rarefied gases at the moment, not for dense configurations.

It made sense to combine two other projects into this one: `fpart` (force-based interactions) and `sphpart` (smoothed particle hydrodynamics). Essentially this amounts to the different methods for integrating the three different Lagrangian systems,

1. `epart`:
$$ L=\sum_i \frac{1}{2}m_i v_i^2-\sum_{i\lt j}\begin{cases} \infty & \text{if $|x_i-x_j|\lt r_{ij}$}\\ 0 & \text{otherwise}\end{cases}$$

2. `fpart` (for short ranged interaction $V$):
$$ L=\sum_i \frac{1}{2}m_i v_i^2-\sum_{i\lt j}V(|x_i-x_j|)$$

3. `sphpart`:
$$L=\sum_i \frac{1}{2}m_i v_i^2-\sum_i m_i\cdot U\left(\sum_j m_j W(x_i-x_j)\right) $$

`fpart` and `sphpart` are not implemented yet.

## epart Notes

Example of 1m particles evolving over time in a free expansion setup:

![1mparticles](git-images/idealgas_preview.gif)

It's "almost" exact because only one collision can be handled per particle per timestep. This is good enough for very rarified gases, but it would be a problem 
for higher densities. Maybe we can drop in an exact event driven molecular dynamics code for the solver later.

