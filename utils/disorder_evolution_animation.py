import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from math import gcd

class DisorderEvolution:
  """Animate the evolution of Hofstadter butterflies with increasing disorder for different lattice types"""

  def __init__(self, lattice_class, length=20, t=1.0, max_W=2.0, q=50,
                 frames=51, interval=100, style='dark_background'):

    """
    Initialize DisorderEvolution class.

    Parameters:
      lattice_class: Class object for the lattice type (Square_Hamiltonian, Triangular_Hamiltonian, etc.)
      length (int): Lattice size (default: 20)
      t (float): Hopping parameter (default: 1.0)
      max_W (float): Maximum disorder strength (default: 2.0)
      q (int): Maximum denominator for phi values (default: 50)
      frames (int): Number of frames for animation (default: 51)
      interval (int): Time between frames in milliseconds (default: 100)
      style (str): Matplotlib style to use (default: 'dark_background')
    """
    self.lattice_class = lattice_class
    self.length = length
    self.t = t
    self.max_W = max_W
    self.q = q
    self.frames = frames
    self.interval = interval

    # Animation properties
    self.fig = None
    self.ax = None
    self.scatter = None
    self.title = None

    # Set plot style
    plt.style.use(style)

    # Initialize system without disorder
    self.model = self.lattice_class(
        length=self.length,
        t=self.t,
        W=0.0,
        phi=0.0,
        q=self.q
    )

  def _collect_butterfly_data(self, W):
    """
    Collect Hofstadter butterfly data for a given disorder strength.

    Parameters:
      W (float): Disorder strength

    Returns:
      tuple: Arrays of phi values and corresponding energies
    """
    model = self.lattice_class(
        length=self.length,
        t=self.t,
        W=W,
        phi=0.0,
        q=self.q
    )

    phis = []
    energies = []

    for q_val in range(1, model.max_q + 1):
      for p in range(q_val + 1):
        if gcd(p, q_val) == 1:
          model.phi = p / q_val
          model.construct_hamiltonian()
          phis.extend([model.phi] * model.N)
          energies.extend(model.evals.tolist())

    return np.array(phis), np.array(energies)

  def _init_plot(self):
    """Initialize the animation plot"""
    self.fig, self.ax = plt.subplots(figsize=(12, 8))
    self.ax.set_xlabel('Flux per Plaquette Φ')
    self.ax.set_ylabel('Energy E')

    # Get initial data
    phis, energies = self._collect_butterfly_data(W=0.0)

    # Create scatter plot
    self.scatter = self.ax.scatter(phis, energies, s=0.1, color='blue', alpha=0.5)
    self.title = self.ax.set_title(
        f'Hofstadter Butterfly for {self.model.lattice_type} Lattice (W = 0.00)'
    )

    # Set plot limits and grid
    self.ax.set_xlim(-0.1, 1.1)
    self.ax.set_ylim(-4, 4)
    self.ax.grid(True, alpha=0.2)

  def _update(self, frame):
    """Update function for animation frame"""
    W = (frame / (self.frames - 1)) * self.max_W

    # Collect new data
    phis, energies = self._collect_butterfly_data(W)

    # Update scatter plot
    self.scatter.set_offsets(np.c_[phis, energies])
    self.title.set_text(
        f'Hofstadter Butterfly for {self.model.lattice_type} Lattice (W = {W:.2f})'
    )

    return self.scatter, self.title

  def animate(self, save_path=None):
    """
    Generate and display/save the disorder evolution animation.

    Parameters:
      save_path (str, optional): Path to save the animation file
    """
    self._init_plot()

    anim = FuncAnimation(
        self.fig,
        self._update,
        frames=self.frames,
        interval=self.interval,
        blit=True
    )

    if save_path:
      print(f"Saving animation to {save_path}")
      anim.save(save_path, writer='pillow', fps=10)
    else:
      plt.show()

    plt.close()

# Example usage:
if __name__ == "__main__":
  # Import any lattice type
  # For example, the square
  from square import Square_Hamiltonian

  # Create animation for square lattice
  evolution = DisorderEvolution(
      lattice_class=Square_Hamiltonian,
      length=20,
      t=1.0,
      max_W=2.0,
      q=50
  )

  # Generate and save animation
  evolution.animate(save_path='square_butterfly_evolution.gif')

# Examples for other lattice types:
# from triangular import Triangular_Hamiltonian
# from kagome import Kagome_Hamiltonian
# from honeycomb import Honeycomb_Hamiltonian
