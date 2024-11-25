# Simulating Anderson Localization and Hofstader Butterflies on Common Euclidean Lattice Structures
 ## Created by Amal Bumbia, Rida Siddiqi, and Max Wroblewski.

# Getting Started

To use the notebooks, ....
- talk about repo structure breifly

# Introduction

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
