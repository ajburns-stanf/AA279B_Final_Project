# Low-ΔV Earth-to-L2 Halo Transfers via Stable Manifolds

AA 279B final project — Alexander Burns, Kaila Coimbra, Yuji Takubo

MATLAB framework for designing low-ΔV transfers from low Earth orbit to an Earth-Moon L2 halo orbit in the circular restricted three-body problem (CR3BP).

## Approach

1. **Halo orbit and stable manifold:** Differentially correct an L2 halo orbit, use the monodromy matrix to find its stable eigendirection, and propagate perturbed states backward to generate the stable manifold.
2. **Earth-to-manifold transfer:** Grid search over halo phase and perturbation direction for manifold arcs that pass closest to Earth, then use differential correction to solve for the impulsive maneuvers connecting a given LEO to the manifold.
3. **Optimization:** Apply primer vector theory to check optimality and adjust maneuver timing to reduce total ΔV.

## References

1. D. P. Gordon, "Transfers to Earth-Moon L2 Halo Orbits," M.S. thesis, Purdue University, 2008.
2. G. Bucchioni et al., "Optimal time-fixed impulsive non-Keplerian orbit to orbit transfer algorithm based on primer vector theory," *Communications in Nonlinear Science and Numerical Simulation* 124 (2023): 107307.
