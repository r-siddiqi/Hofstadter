# Simulating Anderson Localization and Hofstader Butterflies on Common Euclidean Lattice Structures
 ### Created by Amal Bumbia, Rida Siddiqi, and Max Wroblewski.

# Getting Started

To use the notebooks, ....
- talk about repo structure breifly

# Introduction
Condensed matter physics realizes the importance of 2D materials when it comes to the observation of interesting phenomena (ie. optical, magnetic, etc) and topological properties as well as technological applications. 
Different 2D materials have different underlying lattice structures. Some common ones are the square lattice, the honeycomb lattice, triangular/hexagonal lattice, and the kagome lattice. 

These lattices can be described via Hamiltonians encoding how a species could "hop" between adjacent sites (nearest neighbors). A useful way to understand these systems is by diagonalizing the associated Hamiltonians to obtain allowed energies and states, and use these eigenvalues and eigenvectors to compute various parameters. In general, it is useful to plot the eigenvalue spectrum, different eigenvectors, density of states, and the inverse participation ratio. In the case of a magnetic field application perpendicular to the lattice, plotting the Hofstader butterfly is a primary result. Generally, the band structure of these systems is useful to compute, but energy band theory is dependent on the concept of crytal momentum --- a relevant parameter when it comes to describing translational symmetry in the lattice via Bloch's theorem. However, Block's theorem does not hold in cases like hyperbolic lattice connectivity. There is also little support for computing Hofstader butterflies while there exist multiple Python packages referring to band theory calculations. This brings us to what is interesting about our implementation of lattice simulations in this repository:

- We create base cases for the most common, physically realizable lattice types without reliance on crystal momentum to allow for extending this code to cases where Bloch's theorem no longer holds.
- We account for the presence of additional on-site disorder as well as a constant, perpendicular magnetic field. This allows us to simulate Anderson localization and plot Hofstader butterflies with and without disorder. The magnetic field support also sets the stage for further invesitgation into the quantum hall effect.
- The classes creating the hamiltonians for each lattice type can be modified or refactored to support additional parameters or more accurately reflect a particular type of material such as graphene.



# Basic Usage

# Plot Examples

Eigenvalue spectrum:

Eigenvector:

Inverse participation ratio:

Density of states:

Square butterfly with no disorder:
Honeycomb butterfly with no disorder:

An interesting result discussed in [PAPER] is that the presence of disorder kills the butterfly structure.

Square butterfly with marginal disorder:

Square butterfly with strong disorder:

Many other example plots can be found in the "plots" folder.

# Potential Applications
The code in this repository simulates some basic tight-binding hamiltonians specific to common 2D lattice types found in materials with the option to introduce on-site disorder or a uniform magnetic field. Modifications would need to be made in order to simulate specific 2D materials or account for more complicated/specific defects, a non-uniform magnetic field, or other complex systems. 
A non-exhaustive list of possible applications of this code is as follows: 
- Simulating graphene (and possibly twisted-bilayer graphene)
- Hyperbolic Circuit QED

# Resources
Here are some papers we referenced and loosely recreated the results of. You may also find them interesting or useful to read.

# Future Directions
Some possible modifications to this code could involve:
- Support for hyperbolic cases
- Support for open boundary conditions
- Extensions regarding parameters useful to understanding the quantum hall effect on such systems as well as topological properties (ie. hall conductance, thouless conductance, chern number)

# Credits
This started as a project for PHY 381C (Computational Physics) at UT Austin (class website). The presentation we gave in-class on this repository is in the "presentations" folder.
Special thanks to Dr. William Gilpin (wgilpin@utexas.edu) for being a great instructor all semester!
