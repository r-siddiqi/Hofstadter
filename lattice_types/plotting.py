import numpy as np
import matplotlib.pyplot as plt
import os

class Plot_Lattice:
    
    def __init__(self, outputs, save = False):
        
        
        self.L = outputs[0]  # Lattice dimension
        self.N = self.L * self.L  # Total number of sites
        self.t = outputs[1]  # Hopping parameter
        self.disorder = outputs[2]  # Disorder strength
        self.phi = outputs[3]  # Magnetic flux per plaquette
        self.max_q = outputs[4]  # Maximum denominator for phi values
        self.evals = outputs[5]
        self.evecs = outputs[6]
        self.lattice_type = outputs[7]
        self.save = save


        if self.save == True:

            if self.disorder == 0:
                self.path = 'plots\\' + self.lattice_type + '\\No_Disorder'
            else:
                self.path = 'plots\\' + self.lattice_type + '\\Disorder'
            
            if not os.path.exists(self.path):
                os.makedirs(self.path)
    
    
    """ Basic plotting functions """

    def plot_evals(self):
        # Plot eigenvalues of hamiltonian matrix
        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.plot(self.evals, '.')
        plt.ylabel(r'Eigenvalues $E_i$')
        plt.xlabel('Index $i$')
        plt.title('Eigenvalues of the ' + self.lattice_type + ' Lattice Hamiltonian')
        plt.legend([legend])
        plt.grid(True)
        plt.show()

    def plot_evec(self):
        # Plot some eigenvector in the middle of the spectrum
        self.psi = self.evecs[:,self.L//2]

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.plot(np.abs(self.psi)**2)
        plt.xlabel('x')
        plt.ylabel(r'$ |\psi(x)|^2$')
        plt.title('Arbitrary Eigenvector')
        plt.legend([legend])
        plt.grid(True)
        plt.show()

    def plot_evec_disorder(self):
        # Plot some eigenvector in the middle of the spectrum in the presence of disorder
        self.psi = self.evecs[:,self.L//2] # Some eigenvector in the middle of the spectrum

        fig, ax = plt.subplots(2,1,sharex=True)
        ax[0].plot(np.abs(self.psi)**2)
        ax[1].semilogy(np.abs(self.psi)**2)
        ax[1].set_xlabel('x')
        ax[0].set_ylabel(r'$ |\psi(x)|^2$')
        ax[1].set_ylabel(r'$ |\psi(x)|^2$')

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.title('Arbitrary Eigenvector')
        plt.legend([legend])
        plt.grid(True)
        plt.show()

    def plot_pr(self):
        # Plot Participation Ratio
        self.PR = 1./np.sum(np.abs(self.evecs)**4, axis=0) # 'evecs' is a matrix of $\psi_i(x)$ amplitudes, 1st axis is x. This does the sum over x.

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.plot(self.evals, self.PR, 'o')
        plt.xlabel('Energy $E$')
        plt.ylabel('Inverse Participation Ratio (IPR)')
        plt.title('Localization Properties of the '+self.lattice_type+' Lattice')
        plt.legend([legend])
        plt.grid(True)
        plt.show()

    """ Advanced computations and plotting """

    def plot_density_of_states(self, sigma=0.1, num_points=1000):
        """
        Plot the density of states.

        Parameters:
            sigma (float): Standard deviation for Gaussian broadening.
            num_points (int): Number of points in the energy grid.
        """
        energy_min = np.min(self.evals) - 1
        energy_max = np.max(self.evals) + 1
        E_vals = np.linspace(energy_min, energy_max, num_points) # Artificial energy space seperate to eigenvalues
        dos = np.zeros_like(E_vals)

        for E_n in self.evals:
            dos += np.exp(-((E_vals - E_n) ** 2) / (2 * sigma ** 2)) / (np.sqrt(2 * np.pi) * sigma) # Using gaussian broadening

        legend = f'L={self.L}, t={self.t}, W={self.disorder}, $\phi$={self.phi}'
        plt.figure(figsize=(8, 6))
        plt.plot(E_vals, dos)
        plt.xlabel('Energy $E$')
        plt.ylabel('Density of States $g(E)$')
        plt.title('Density of States vs Energy for '+ self.lattice_type + ' Lattice')
        plt.legend([legend])
        plt.grid(True)
        plt.show()
