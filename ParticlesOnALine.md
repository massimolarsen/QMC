# Particles on a Line
To test out the VEGAS algorithm we will first apply the algorithm to a simplified transport problem.  The things to do are:

1. [ ] Implement a standard Monte Carlo solution to the transport problem as described below to solve for the leakage out of the left and right sides of the slab.
2. [ ] Apply the VEGAS algorithm using the sampled source location as a biasing parameter 
    - [ ] Optimize/solve for the leakage to the left and right side of the slab independently
    - [ ] Optimize/solve for the total leakage from the slab
    - [ ] Optimize/solve for the leakage to the left and right side of the slab and total leakage in a single integral

## Problem Description
Consider a simplified one-speed transport problem where there is only one spatial dimension $x$ with a direction cosine $\mu$, and particles can either go forward $\mu = 1$ or backward $\mu = -1$.  Define the angular flux forward and backward as $\psi^+$ and $\psi^-$, respectively (in units of per cosine cm$^{-2}$).  The slab where we are interested in the solution is from $x=[0,L]$, i.e., it has a width $L$ (cm).  Within the slab, the macroscopic absorption cross section is $\sigma_a$ (cm$^{-1}$), the scattering cross section is $\sigma_s$, and the total cross section is $\sigma_t = \sigma_a + \sigma_s$.  The scattering probability is isotropic.  There is an extraneous, isotropic source $q(x)$ (particles-cm$^{-3}$) that is a constant value $q_o$ over the left $W$ cm of the slab, with $0 < W \leq  L$.  There is vacuum boundary conditions, i.e., $\psi^+(0)=0$ and $\psi^-(L)=0$.

## Monte Carlo Process.
The MC analog transport process should be implemented with the following algorithm for $N$ independent histories:

1. Sample a source particle uniformly between 0 and $W$, with direction $\mu=1$ with probability 1/2 and $\mu=-1$ with probability 1/2.
2. Sample a distance to collision $x_c$, from $p(x_c) = \sigma_t \exp (-\sigma_t x_c), \quad x_c > 0$
3. If $|x_c|$ is greater than the distance to a boundary (depending on direction), then the particle has leaked.  Add the weight (1.0 for analog) to the correct leakage tally for the right (+) or left (-) boundary, terminate the history, and return to the beginning of the algorithm.  Since there are vacuum boundary conditions and the particle only has unit direction, the surface current tallies (i.e., estimating the current $J^+(L)$ and $J^-(0)$) are simply the mean of the weights of particles that leak, e.g.,
$$
   J^+ = \frac{1}{N} \sum_{i=1}^M  w(L) 
$$
where $M$ is the number of particles that crossed the right boundary with weight $w(L)$ (for analog this is just 1.0, but we will need to adjust this for the VEGAS algorithm possibly).
4. If $x_c$ is less than the distance to a boundary, the particle has collided.  Move the particle to the new position as $x += \mu x_c$, where $\mu=\pm 1$.
5. Sample a random number $\eta \sim Unif(0,1)$. 
	1. If $\eta < \sigma_a/\sigma_t$, then the particle is absorbed.  Terminate the history and return to beginning of algorithm
	2. Else: The particle has scattered
		1. Sample *another* random number $\eta$.  If $\eta < 0.5$, the particle has backscattered, so change the direction from $\mu = \pm1$ to $\mu = \mp 1$.
		2. Else, the particle "forward" scatters, so it doesn't do anything.  Return to step $2$ of the algorithm  

To get the current in the correct units, you can multiply the tally results by $q_0$.
The forward scattering process is a bit unusual, but it is a result of how we defined the scattering cross section.  You could have instead defined only the backscattering portion, but it will make the equations unusual compared to normal transport equations.

## VEGAS Algorithm
TODO

## Reference Solution
There are analytic solutions to this problem, but only for uniform source throughout, which doesn't give any interesting shape to the problem to challenge the VEGAS algorithm, or for no source with boundary conditions, which has nothing for the VEGAS algorithm to optimize.  Instead, we will derive a deterministic numerical solution for the problem described above using a finite-difference method. 

Since we are in units of per cosine, we can define the scalar flux as $\phi(x) = \psi^+(x) + \psi^-(x)$.  Simplifying the Boltzmann transport equation to account for the one spatial direction and two directions, the governing equations become only a function of $x$ as
$$
\frac{d \psi^+(x)}{d x} + \sigma_t \psi^+ = \frac{\sigma_s}{2}\left( \psi^+(x) + \psi^-(x) \right) + \frac{q(x)}{2}
$$
and
$$
-\frac{d \psi^-(x)}{d x} + \sigma_t \psi^- = \frac{\sigma_s}{2}\left( \psi^+(x) + \psi^-(x) \right) + \frac{q(x)}{2}.
$$

To apply the finite difference method in the spatial dimension (no need for the angular dimension here), we will integrate the equations over a cell of width $\Delta x$ with edges $x_{i-1/2}$ and $x_{i+1/2}$, for $i=0,\ldots,K$.  To close the equations, we will apply the diamond-difference discretization, where we assume the flux, in each direction, is continuous and the average of the flux over the cell is just the linear average of the two edges, i.e., $\psi_i = (\psi_{i-1/2} + \psi_{i+1/2})/2$. The integration of the governing equations with the linear approximation gives, for the forward direction, the following equation:
$$
\psi^+_{i+1/2} - \psi^+_{i-1/2} + \frac{\sigma_t \Delta x}{2}\left(\psi^+_{i+1/2} + \psi^+_{i-1/2}\right) - \frac{\sigma_s \Delta x}{4}\left[\psi^+_{i+1/2} + \psi^+_{i-1/2} + \psi^-_{i+1/2} + \psi^-_{i-1/2} \right] = \frac{q_0 \Delta x}{2}
$$
for cells where $q(x)$ is non-zero, otherwise there is no $q_0$ term.  The equations for the negative direction are
$$
-\psi^-_{i+1/2} + \psi^-_{i-1/2} + \frac{\sigma_t \Delta x}{2}\left(\psi^-_{i+1/2} + \psi^-_{i-1/2}\right) - \frac{\sigma_s \Delta x}{4}\left[\psi^+_{i+1/2} + \psi^+_{i-1/2} + \psi^-_{i+1/2} + \psi^-_{i-1/2} \right] = \frac{q_0 \Delta x}{2}
$$


The boundary conditions are applied (I think) by defining for the positive direction that the incident current is 0, which is a bit confusing because $\mu=1$, but simplifies as follows:
$$
\int_{x_{-1/2}}^{x_{+1/2}} dx \; J^+(x) = \int_0^1 d \mu \int_{x_{-1/2}}^{x_{+1/2}} dx \frac{d \psi^+}{d x}\delta(\mu-1) = J^+_{1/2} - J^+_{-1/2} =  \psi^+_{1/2} - 0
$$
Applying the same integration and approximations as before then gives the following equation for cell $i=0$:
$$
\psi^+_{1/2} + \frac{\sigma_t \Delta x}{2}\left(\psi^+_{1/2} + \psi^+_{-1/2}\right) - \frac{\sigma_s \Delta x}{4}\left[\psi^+_{1/2} + \psi^+_{-1/2} + \psi^-_{1/2} + \psi^-_{-1/2} \right] = \frac{q_0 \Delta x}{2}
$$
There is a similar expression derived by $J^-_{K+1/2}=0$ for the equation in the negative direction of the last cell.  (It is possible I will need to also make each of the $\psi^+_{-1/2}$ terms zero in the above equation since $J^+ = \psi^+$, but I will have to play around with it a little to make sure I did this part correctly).

This should give us $2K+4$ equations for $2K+4$ edge unknowns ($K+2$ unknowns in each direction because we have $K+1$ cells), which we can assemble into a global matrix and solve directly.  We can then compare the edge angular fluxes to the currents solved for with Monte Carlo.  

 