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

# Potential Applications

# Resources

# Future Directions

# Credits

In condensed matter research, you often need to run simulations. There is a notable lack of open source code for base cases that you can modify to suit your own purposes. There are Python packages that have attempted to solve this issue, but they are difficult to customise and use. 

This code is cleanly compartmentalized into classes featuring common lattice structures, many of which are found in 2D materials with promising properties under current study. It is much easier to modify a class of a specific lattice type and compute only what you need.

objectives:
- make the final plots nicer, more clear 
- complete cases for kagome, oblique, triangular
- find a clean way to account for redundancy between classes - for instance, many of the plotting functions wind up being the same for some (if not all) of the classes. We might want to make them global and call the plotting outside of the class
- make a folder of all example runs and plots
- add in error statements in case of bad user inputs
- fix comments and docstrings

interesting observations:
- as you play around with changing the initial parameters, notice that adding disorder to the system kills the butterfly structure
