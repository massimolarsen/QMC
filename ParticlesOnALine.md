# Particles on a Line
To test out the VEGAS algorithm we will first apply the algorithm to a simplified transport problem.  The things to do are:

1. [ ] Implement a standard Monte Carlo solution to the transport problem as described below to solve for the leakage out of the left and right sides of the slab.
2. [ ] Apply the VEGAS algorithm using the sampled source location as a biasing parameter 
    - [ ] Optimize/solve for the leakage to the left and right side of the slab independently
    - [ ] Optimize/solve for the total leakage from the slab
    - [ ] Optimize/solve for the leakage to the left and right side of the slab and total leakage in a single integral

## Problem Description
Consider a simplified one-speed transport problem where there is only one spatial dimension $x$ with a direction cosine $\mu$, and particles can be constrained to scatter isotropically or further simplified to scatter only forward $\mu = 1$ and backward $\mu = -1$.  Define the angular flux forward and backward as $\psi^+$ and $\psi^-$, respectively (in units of per cosine cm$^{-2}$).  The slab where we are interested in the solution is from $x=[0,L]$, i.e., it has a width $L$ (cm).  Within the slab, the macroscopic absorption cross section is $\sigma_a$ (cm$^{-1}$), 

the scattering cross section is $\sigma_s$, and the total cross section is $\sigma_t = \sigma_a + \sigma_s$. There is an extraneous, isotropic source $q(x)$ (particles-cm$^{-3}$) that is a constant value $q_o$ over the left $W$ cm of the slab, with $0 < W \leq  L$.  There is vacuum boundary conditions, i.e., $\psi^+(0)=0$ and $\psi^-(L)=0$.

## Monte Carlo Process.
The MC analog transport process should be implemented with the following algorithm for $N$ independent histories:

1. Sample a source particle uniformly between 0 and $W$, with an isotropically sampled angle $\mu$ 
    1. In the simple problem we further constrain the particles angle to be $\mu=1$ with probability 1/2 and $\mu=-1$ with probability 1/2.
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

To get the current in the correct units, we multiply the tally results by $q_0$.

## VEGAS Algorithm


### Importance Sampling
The VEGAS algorithm is based on a standard Monte Carlo biasing technique called importance sampling, which works by 
adjusting the distribution from which the Monte Carlo function samples from. Ideally it will cause the Monte Carlo function to sample more often from regions of more "importance". Given an integral $I$ and a function $f(x)$ over some domain $D$ we can construct the importance sampling integral by simply multiplying by 1, where 1 = $\frac{f^*(x)}{f^*(x)}$.

$$
I = \int_D \frac{f(x)}{f^*(x)} f^*(x) dx
$$

Given that $f^*(x) \neq 0$ this leads to the Monte Carlo estimate:

$$
\tilde{f}(X) = \frac{1}{N} \sum_{i=1}^N \frac{f(X_i)}{f^*(X_i)}
$$

The variance of $f(X)$ can be shown to be minimized when $f^*(x) \propto |f(x)|$ and if $f^*(x) = \frac{|f(x)|}{\int_D |f(x)|}$ the variance disappears and the error is 0. This exact distribution is not pratical, or even possible, to sample from in reality so importance sampling algorithms, such as VEGAS, attempt to approximate this distribution as efficiently and accurately as possible.

### VEGAS

VEGAS's main goal is to minimize the statistical error of the Monte Carlo estimate of $I$. However without knowing the analytical solution, the exact $f^*$ can not be calculated in one step or even at all. In order to accurately estimate $f^*$ VEGAS employs an iterative technique where it remaps the integration variables of $f^*$ every iteration in an attempt to minimize the error.

Given an integration variable $x$ VEGAS would attempt to transform it into a new variable $y$ using a grid created in the original $x$ space. The variance for the Monte Carlo estimate of this transformation is:

$$
Var(I) = \frac{1}{M} (\sum_i J_i \int_{x_i}^{x_{i+1}} dx f^2(x) - I^2)
$$

Where $J_i$ is the Jacobian of the transformation function.

$$
J_i = J(y) = N \Delta x_i
$$

By constraining the Jacobian's to follow:

$$
\sum_i \frac{\Delta x_i}{J_i} = \sum_i \Delta y_i = 1
$$


###### (need to derive)
It can be shown that :

$$
\frac{J^2_i}{\Delta x_i} \int_{x_i}^{x_{i+1}} dx f^2(x) = N^2 \Delta x_i \int_{x_i}^{x_{i+1}}dx f^2(x) = constant
$$

VEGAS uses this final condition as the basis for its iterative remapping by adjusting the grid until the above condition is satisfied.

After creating the grid VEGAS approximates a new $f^*$ with a set of histograms to create a sampling distribution for the following pass. VEGAS also ensures that each newly created bin has an equiprobable chance of being sampled. Such a tactic reduces the memory needs of the algorithm. Once the per-iteration error stops decreasing, VEGAS has acheived an accurate mapping of the $f^*$ function and can then produce a reasonable Monte Carlo estimate of the initial integrand $f$. 

### Preconditioning
Since VEGAS begins each set of iterations with no information on the integrand it is common for the initial set of estimates to be wildy incorrect. In that case it is possible to first "train" the algorithm on a throwaway set of iterations, discard the results, and then run the real evalution. This method can even be adapted to more complex cases where the algorithm can be trained on a simpler radiography example and then applied to the more complex problem. This technique will reduce run time as well as produce more accurate estiamtes.

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


## Tests
### Spatial Test

We ran two practice problems to test the abilities of the VEGAS algorithm and learn more about how the VEGAS python package works. In our test function we allowed VEGAS to control both the spatial position ($x$) as well as the angular component ($\mu$) of the source particles. After the particles were spawned in by VEGAS they went through a normal Monte Carlo history before returning their tally back to the algorithm. In this case where we were solely looking at total leakage outside of the slab we simply returned to VEGAS a 1 if the particle leaked or a 0 if it was absorbed. Additional tallies, such as left/right edge leakage, were kept but these were not given to VEGAS. This meant that the algorithm was adjusting where it spawned particles and what direction they were heading in an effort to optimize the overall leakage of the slab.

The first test problem used a simple setup where the slab ($L=5$) had an evenly distributed isotropic source ($q_0$) with a strength of 3 neutrons/cm\^3/second. The particles were then further constrained to only the forward and backwards directions where $\mu=1$ and $\mu=-1$. The scattering cross section ($\sigma_s$) and absorption cross sections ($\sigma_a$) were 5b and 0.3b respectively. 

The results from a standard Monte Carlo estimate, the reference solution, and the VEGAS estimate were all compared to ensure that our functions were running correctly. The Monte Carlo estimate was run for 100,000 histories while the VEGAS algorithm ran 10 iterations of 10,000 histories. This number of histories does not achieve extremely accurate results but for our purpose of testing it was plenty. 

After verifying that all of our functions agreed with each other we plotted the VEGAS map to give us a visual representation of what the algorithm is actually accomplishing. Figure 1 shows the VEGAS map for our initial test problem. The x-axis represents the x-axis in the problem while the y-axis is showing us the angular distribution ($\mu$). Every bin shown in the image has and equal chance of being sampled from. This means that VEGAS simply has to adjust the grid, as discussed above, rather than storing and adjusting sampling distributions. VEGAS naturally populates the more important areas of the problem with more bins in an effort to sample from there more often. In this case where we have a relatively absorbing slab the edges of the problem contribute more towards the overall leakage, so VEGAS adjusts the grid to spawn less particles in the center and more on the edges. 

### Angular Test

The second test problem was chosen to focus more on the angular component rather than the spatial. In order to force more attention towards the initial angle the absorption cross section was increased to 3b and the source length was reduced to 3cm. Increasing the absorption rate greatly increases the importance of the initial angle since a particle heading away from the edges was far more likely to be absorbed. Then decreasing the source length almost ensured that no particles would be able to leak out the right side of the slab. Ideally this would make particles spawn on the left edge and head in the $\mu = -1$ direction.


This test was run once with full isotropic scattering and once with the constrained unit scattering as we used above. Figure 2 shows a noticeable difference between the techniques as we would expect. The unit direction scattering forces all $\mu < 0$ to head left and all $\mu > 0$ to head right which creates the distinct zones that we can see in Figure 2a. The isotropic scattering shown in 2b has a much more gradual change in the distribution as $\mu \rightarrow 1$.