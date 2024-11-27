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

## Tight-Binding Hamiltonian
Suppose we take an atom. There are orbitals associated with it that describe regions where electrons could exist -- these are eigenfunctions of the Hamiltonian describing our atom. But what if we configure our atom with other atoms in a crystalline configuration? 
How do we describe the energy of a system where particles are placed in a crystalline configuration?
The following tight-binding Hamiltonian describes lattices in the absence of external interactions while accounting for hopping:
$$H = \omega_0 \sum_{i} a^{\dagger}{i} a{i} - t \sum_{\langle i,j \rangle} (a^{\dagger}{i} a{j} + a^{\dagger}{j} a{i})$$
where $\omega_0$ is the on-site energy. The first sum runs over all lattice sites. The second sum describes hopping between nearest neighbors with a hopping amplitude $t$. It also encodes lattice geometry. 

When we're simulating specific lattice types, it's actually easier to construct a Hamiltonian in terms of a matrix populated by considering the lattice geometry. By this, we mean populating an N x N matrix based on the nearest neighbor hopping as it would occur on a particular lattice.

## Anderson Localization
The key idea behind Anderson Localization is that certain materials can undergo a phase transition from conductor to insulator if the system passes a disorder threshold. Thus, in systems of sufficiently large disorder (such as defected semiconductors) the electronic wavefunction associated with a now spatially localized state becomes localized. This localization influences the aforementioned phase transition. \textbf{In other words, spatial localization of the electronic wavefunction causes a change in the conductance of a highly disordered material.}

The Anderson tight-binding model allows us to describe this phenomenon more effectively. Here, electrons can tunnel between neighboring lattice sites until high disorder in the lattice causes quantum amplitudes associated with the tunneling paths to cancel amongst each other. A localized wavefunction follows. Equivalently, we can say that the incoming wave is scattered by potentials arising from the disorder, and the scattered wavelets destructively interfere as they propel forward at high disorder. This causes an exponential decay of the wavefunction. In this way, \textbf{Anderson localization} can be thought of as an interference phenomenon. 

Experimentally, electron localization has mostly been observed in a 1D case. 

The Anderson Hamiltonian can help us describe the localization in more technical terms. We write it as such

$$H = W \sum_{n} (\varepsilon_{n} c^{\dagger}{n} c{n}) + t\sum_{\langle n,m \rangle} (c^{\dagger}{n} c{m} + \text{h.c.})$$ where $t$ is the parameter describing the nearest hopping neighbor, $W$ is the disorder parameter, and $\epsilon_n$ is the random on-site energy in the range $[-\frac{1}{2},\frac{1}{2}]$.

## Inverse Participation Ratio (IPR)

The "participation ratio" gives an estimation of the localization length:
$$IPR = \frac{(\sum_x |\psi(x)|^2)^2}{ \sum_x |\psi(x)|^4}$$
(the numerator is not necessary if wavefunctions are normalized).

## Hofstader Butterflies
What happens when we apply a perpendicular, uniform magnetic field onto a lattice? The general tight-binding hamiltonian will now involve a "Peierls phase" accounting for the magnetic flux through each plaquette as well as relevant changes in the boundary conditions.
An interesting result is that if we plot the energies as a function of magnetic flux ratios ($\phi = p/q$) such that $p$ and $q$ are coprime integers, we obtain a fractal pattern. It is a recursive structure. The way we constructed the butterfly involved choosing a maximum value for $q$, iterating through all the coprime $p$,$q$ pairs leading up to that point, and then reconstructing the hamiltonian for each consequent $\phi = p/q$. 

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
- Simulating graphene (twisted-bilayer graphene)
- Hyperbolic Circuit QED
- Realizing the quantum hall effect in various moiré materials

# Resources
Here are some papers we referenced or loosely recreated the results of. You may also find them interesting or useful to read.

- https://www.researchgate.net/publication/288374712_Energy_spectrum_of_a_honeycomb_lattice_under_nonuniform_magnetic_fields
- https://physics.bgu.ac.il/ARCHIVE/public_projects/2018_06_25_idofried.pdf
- https://physics.bu.edu/~okctsui/PY543/5_notes_Tight%20Binding.pdf
- https://iopscience.iop.org/article/10.1088/1367-2630/ac4126/pdf
- https://www.researchgate.net/publication/286357286_Energy_Spectrum_of_a_Triangular_Lattice_in_a_Uniform_Magnetic_Field_Effect_of_Next-Nearest-Neighbor_Hopping
- https://journals.jps.jp/doi/10.1143/JPSJ.53.3101#:~:text=A%20computer%20simulation%20is%20performed%20to%20study%20the,Thouless%20number%20g%20%28%20L%20%29%20is%20employed.
- https://link.springer.com/chapter/10.1007/978-3-030-21511-8_5
- https://journals.aps.org/prb/abstract/10.1103/PhysRevB.71.125310
- https://www.degruyter.com/document/doi/10.1515/9781400846733/html
- https://www.sciencedirect.com/science/article/abs/pii/0038109876912539
- https://phas.ubc.ca/~berciu/TEACHING/PHYS502/PROJECTS/17AL.pdf
- http://arxiv.org/abs/2310.07978
- https://www.damtp.cam.ac.uk/user/tong/qhe/qhe.pdf
  

# Future Directions
Some possible modifications to this code could involve:
- Support for hyperbolic cases
- Support for open boundary conditions
- Extensions regarding parameters useful to understanding the quantum hall effect on such systems as well as topological properties (ie. hall conductance, thouless conductance, chern number)

# Credits
This started as a project for PHY 381C (Computational Physics) at UT Austin (class website: https://www.wgilpin.com/cphy/). The presentation we gave in-class on this repository is in the "presentations" folder.
Special thanks to Dr. William Gilpin (wgilpin@utexas.edu) for being a great instructor all semester!
