# Simulating Anderson Localization and Hofstadter Butterflies on Common Euclidean Lattice Structures
 ### Created by Amal Bumbia, Rida Siddiqi, and Max Wroblewski.

# Getting Started

To use the notebooks, ensure you have a `conda` environment with:
- `numpy`
- `matplotlib`

# Introduction
Condensed matter physics realizes the importance of 2D materials when it comes to the observation of interesting phenomena (ie. optical, magnetic, etc) and topological properties as well as technological applications. 
Different 2D materials have different underlying lattice structures. Some common ones are the square lattice, the honeycomb lattice, triangular/hexagonal lattice, and the kagome lattice. 

These lattices can be described via Hamiltonians encoding how a species could "hop" between adjacent sites (nearest neighbors). A useful way to understand these systems is by diagonalizing the associated Hamiltonians to obtain allowed energies and states, and use these eigenvalues and eigenvectors to compute various parameters. In general, it is useful to plot the eigenvalue spectrum, different eigenvectors, density of states, and the inverse participation ratio. In the case of a magnetic field application perpendicular to the lattice, plotting the Hofstader butterfly is a primary result. Generally, the band structure of these systems is useful to compute, but energy band theory is dependent on the concept of crytal momentum --- a relevant parameter when it comes to describing translational symmetry in the lattice via Bloch's theorem. However, Block's theorem does not hold in cases like hyperbolic lattice connectivity. There is also little support for computing Hofstadter butterflies while there exist multiple Python packages referring to band theory calculations. This brings us to what is interesting about our implementation of lattice simulations in this repository:

- We create base cases for the most common, physically realizable lattice types without reliance on crystal momentum to allow for extending this code to cases where Bloch's theorem no longer holds.
- We account for the presence of additional on-site disorder as well as a constant, perpendicular magnetic field. This allows us to simulate Anderson localization and plot Hofstader butterflies with and without disorder. The magnetic field support also sets the stage for further invesitgation into the quantum hall effect.
- The classes creating the hamiltonians for each lattice type can be modified or refactored to support additional parameters or more accurately reflect a particular type of material such as graphene.

## Tight-Binding Hamiltonian
Suppose we take an atom. There are orbitals associated with it that describe regions where electrons could exist -- these are eigenfunctions of the Hamiltonian describing our atom. But what if we configure our atom with other atoms in a crystalline configuration? 
How do we describe the energy of a system where particles are placed in a crystalline configuration?
The following tight-binding Hamiltonian describes lattices in the absence of external interactions while accounting for hopping:

$$H = \omega_0 \sum_i a_i^{\dagger} a_i - t \sum_{\langle i,j\rangle} (a_i^{\dagger} a_j + a_j^{\dagger} a_i)$$

where $\omega_0$ is the on-site energy. The first sum runs over all lattice sites. The second sum describes hopping between nearest neighbors with a hopping amplitude $t$. It also encodes lattice geometry. 

When we're simulating specific lattice types, it's actually easier to construct a Hamiltonian in terms of a matrix populated by considering the lattice geometry. By this, we mean populating an $N$ x $N$ matrix based on the nearest neighbor hopping as it would occur on a particular lattice.

## Anderson Localization
The key idea behind Anderson Localization is that certain materials can undergo a phase transition from conductor to insulator if the system passes a disorder threshold. Thus, in systems of sufficiently large disorder (such as defected semiconductors) the electronic wavefunction associated with a now spatially localized state becomes localized. This localization influences the aforementioned phase transition. \textbf{In other words, spatial localization of the electronic wavefunction causes a change in the conductance of a highly disordered material.}

The Anderson tight-binding model allows us to describe this phenomenon more effectively. Here, electrons can tunnel between neighboring lattice sites until high disorder in the lattice causes quantum amplitudes associated with the tunneling paths to cancel amongst each other. A localized wavefunction follows. Equivalently, we can say that the incoming wave is scattered by potentials arising from the disorder, and the scattered wavelets destructively interfere as they propel forward at high disorder. This causes an exponential decay of the wavefunction. In this way, \textbf{Anderson localization} can be thought of as an interference phenomenon. 

Experimentally, electron localization has mostly been observed in a 1D case. 

The Anderson Hamiltonian can help us describe the localization in more technical terms. We write it as such

$$H = W \sum_n (\epsilon_n c_n^\dagger c_n) + t \sum_{<n,m>} (c_n^\dagger c_m + h.c)$$

where $t$ is the parameter describing the nearest hopping neighbor, $W$ is the disorder parameter, and $\epsilon_n$ is the random on-site energy in the range $[-\frac{1}{2},\frac{1}{2}]$.

## Inverse Participation Ratio (IPR)

The "participation ratio" gives an estimation of the localization length:
$$IPR = \frac{(\sum_x |\psi(x)|^2)^2}{ \sum_x |\psi(x)|^4}$$
(the numerator is not necessary if wavefunctions are normalized).

## Hofstadter Butterflies
What happens when we apply a perpendicular, uniform magnetic field onto a lattice? The general tight-binding hamiltonian will now involve a "Peierls phase" accounting for the magnetic flux through each plaquette as well as relevant changes in the boundary conditions.
An interesting result is that if we plot the energies as a function of magnetic flux ratios ($\phi = p/q$) such that $p$ and $q$ are coprime integers, we obtain a fractal pattern. It is a recursive structure. The way we constructed the butterfly involved choosing a maximum value for $q$, iterating through all the coprime $p$, $q$ pairs leading up to that point, and then reconstructing the hamiltonian for each consequent $\phi = p/q$. 

_Hoftsadter butterfly for square lattice as q increases (i.e. magnetic flux decreases)_  
<img src="https://github.com/user-attachments/assets/14c65aef-e5b3-47cf-a517-78489e49d2da" width="525" alt="butterfly_q_evolution">


_Hofstadter butterfly for square lattice as hopping parameter, t, increases_  
<img src="https://github.com/user-attachments/assets/f5e795b0-c945-4492-9c41-7f85e0763d1a" width ="525" alt="butterfly_t_evolution">


# Basic Usage

# Plot Examples

**Eigenvalue spectrum:**

For Square Lattice:

<img width="465" alt="Eig Value Spectrum" src="https://github.com/user-attachments/assets/0269271a-cde6-4fc8-81e3-6ac03dad5242">


**Eigenvector:**

For Square Lattice:

<img width="465" alt="Arbitrary Eig" src="https://github.com/user-attachments/assets/78a9d79f-ea1c-4f89-a17d-871424ecd523">



**Inverse participation ratio:**

For Square Lattice:

<img width="465" alt="IPR" src="https://github.com/user-attachments/assets/cfd7c401-e087-441f-a0f7-d825b4b29368">



**Density of states:**

For Square Lattice:

<img width="465" alt="density_of_states" src="https://github.com/user-attachments/assets/6587e584-cb43-4ee6-b36e-20bd98e49c1a">


**Hofstadter Butterfly:**

Square butterfly with no disorder:

<img width="465" alt="square_no_disorder" src="https://github.com/user-attachments/assets/9b55e5c7-d2b2-49eb-8d85-370bb2e4e94d">


Honeycomb butterfly with no disorder:

<img width="465" alt="honeycomb_no_disorder" src="https://github.com/user-attachments/assets/d7ef1d1c-087b-4b6c-ac2f-0d8ce98c7784">



An interesting result discussed in the [literature](https://link.springer.com/article/10.1140/epjb/e2016-70593-4) is that the presence of disorder kills the butterfly structure (but in the high-disorder limit, some butterfly-like structure may still persist).


Square butterfly with increasing disorder (butterfly-like pattern persists up to ~ _W_ = 1.40):

<img src="https://github.com/user-attachments/assets/5bb791fc-4374-472c-bc28-10767e56497c" width="525" alt="disorder_evolution">

_Disorder effects on other common latice types:_  
<img src="https://github.com/user-attachments/assets/74bbc898-361f-4459-8768-3aba7c7ff649" width="525" alt="triangular_disorder_evolution">  

<img src="https://github.com/user-attachments/assets/0a6f8fcc-b9bc-47fb-8167-a84435d560f4" width="525" alt="kagome_butterfly_evolution">


Many other example plots can be found in the ["plots"](https://github.com/r-siddiqi/Hofstadter/tree/main/plots) folder.

# Potential Applications
The code in this repository simulates some basic tight-binding hamiltonians specific to common 2D lattice types found in materials with the option to introduce on-site disorder or a uniform magnetic field. Modifications would need to be made in order to simulate specific 2D materials or account for more complicated/specific defects, a non-uniform magnetic field, or other complex systems. 
A non-exhaustive list of possible applications of this code is as follows: 
- Simulating graphene (twisted-bilayer graphene)
- Hyperbolic Circuit QED
- Realizing the quantum hall effect in various moiré materials

# Resources
Here are some papers we referenced or loosely recreated the results of. You may also find them interesting or useful to read.

## General Theory
- [Tight Binding Method and the Electronic Energy Band Structure](https://physics.bu.edu/~okctsui/PY543/5_notes_Tight%20Binding.pdf)
- [The Quantum Hall Effect](https://www.damtp.cam.ac.uk/user/tong/qhe/qhe.pdf)  
  David Tong's comprehensive notes on QHE
- [Anderson Localization and Its Ramifications](https://www.degruyter.com/document/doi/10.1515/9781400846733/html)  
  Fundamental theory and applications
- [Low-field electron localization in a magnetic field](https://www.sciencedirect.com/science/article/abs/pii/0038109876912539)  
  Kaveh & Mott's seminal paper on electron localization

## Square Lattice
- [Universal properties of the two-dimensional Anderson model in a magnetic field](https://journals.jps.jp/doi/10.1143/JPSJ.53.3101)  
  Computer simulation study of Anderson localization with detailed analysis of scaling properties
- [Quantum localization on the square lattice](https://phas.ubc.ca/~berciu/TEACHING/PHYS502/PROJECTS/17AL.pdf)  
  Comprehensive analysis of localization properties including numerical methods

## Honeycomb Lattice
- [Energy spectrum of a honeycomb lattice under nonuniform magnetic fields](https://www.researchgate.net/publication/288374712_Energy_spectrum_of_a_honeycomb_lattice_under_nonuniform_magnetic_fields)  
  Analysis of nonuniform field effects on honeycomb lattices

## Triangular/Hexagonal Lattice
- [Energy Spectrum of a Triangular Lattice in a Uniform Magnetic Field](https://www.researchgate.net/publication/286357286_Energy_Spectrum_of_a_Triangular_Lattice_in_a_Uniform_Magnetic_Field_Effect_of_Next-Nearest-Neighbor_Hopping)  
  Effects of next-nearest-neighbor hopping and magnetic field interactions

## Kagome Lattice
- [Floquet Hofstadter butterfly on the kagome and triangular lattices](https://journals.aps.org/prb/pdf/10.1103/PhysRevB.98.245145)  
  Figure 2 was used to verify our results for the kagome lattice
- [Hofstadter Butterfly and Many-Body Effects in the Kagome Lattice](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.71.125310)  
  Detailed analysis of magnetic field effects and band structure
- [Electronic Structure of the Kagome Lattice](https://iopscience.iop.org/article/10.1088/1367-2630/ac4126/pdf)  
  Comprehensive band structure analysis and topological properties
  

# Future Directions
Some possible modifications to this code could involve:
- Support for hyperbolic cases
- Support for open boundary conditions
- Extensions regarding parameters useful to understanding the quantum hall effect on such systems as well as topological properties (ie. hall conductance, thouless conductance, chern number)

# Credits
This started as a project for PHY 329 (Computational Physics) at UT Austin ([class website](https://www.wgilpin.com/cphy/)). The presentation we gave in-class on this repository is in the "presentations" folder.
Special thanks to Dr. William Gilpin (wgilpin@utexas.edu) for being a great instructor all semester!
