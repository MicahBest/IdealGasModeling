# IdealGasModeling
### Objective
- Create a general property calculator model for an ideal gas mixture of arbitrary composition of specific components.

$$ (O_2, N_2, CO_2, CO, H_2O, C_x H_y, etc.) $$

- Compute $P, v, T, u, h, s, k, \mu$
- Allow for evaporation/condensation of $H_2O$

### Model Structure
- Table/Matrix that contains the parameters and calculations of the mixture
- Fixed vector length related to the possible species in a mixture (finite list)
- Fixed field width related to the computations needed to compute thermodynamic properties
- Inputs:
  - $P, T$, composition
- Outputs:
  - $M, R, v, u, h, s, k, \mu$
- Model Parameters:
-   $R_u, M, \frac{\bar{C_p}}{R} = f(T), h_f^o$
