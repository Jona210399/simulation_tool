A wrapper tool around pySIMsalabim for jV Curve and EQE simulations.


## Randomization Insights:
- randomizing R_shunt and R_series results in a major decrease of simulation convergence rate. Not randomizing them seems to be the better option since the can be in corporated into the calculation afterwards.

- adding interface traps (setting N_int > 0 ) results in lower convergence rate, but has to be considered. Since the simulated interface trapping depends on the product of C_n/p and N_t_int. So only one of them needs to be randomized. The other one can be set to 1.0. Interface trap energies (E_t_int) need to be in between the energy levels of the interfacing layers.

- removing the layer factor to calculate the charge carrier mobilities (mu_n/p) also results in a worse convergence rate.