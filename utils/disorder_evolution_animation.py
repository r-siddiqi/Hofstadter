import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
plt.style.use('dark_background')

def animate_butterfly_evolution(length=20, t=1.0, max_W=2.0, q=50, save_path=None):
  """
  Animate the evolution of the Hofstadter butterfly as disorder increases

  Args:
      length (int): Lattice size (default: 20)
      t (float): Hopping parameter (default: 1.0)
      max_W (float): Maximum disorder strength (default: 2.0)
      q (int): Maximum denominator for phi values (default: 50)
      save_path (str, optional): Path to save animation
  """
  fig, ax = plt.subplots(figsize=(12, 8))
  ax.set_xlabel('Flux per Plaquette Φ')
  ax.set_ylabel('Energy E')

  # Initial setup with W = 0
  model = Square_Hamiltonian(length=length, t=t, W=0.0, phi=0.0, q=q)
  model.construct_hamiltonian()

  # Initialize plot
  phis = []
  energies = []

  # This follows our plot_hofstadter_butterfly() implementation
  for q_val in range(1, model.max_q + 1):
    for p in range(q_val + 1):
      if gcd(p, q_val) == 1:
        model.phi = p / q_val
        model.construct_hamiltonian()
        phis.extend([model.phi] * model.N)
        energies.extend(model.evals.tolist())

  scatter = ax.scatter(phis, energies, s=0.1, color='blue', alpha=0.5)
  title = ax.set_title('Hofstadter Butterfly for Square Lattice (W = 0.00)')

  ax.set_xlim(-0.1, 1.1)
  ax.set_ylim(-4, 4)
  ax.grid(True, alpha=0.2)

  def update(frame):
    # Convert frame to disorder value (0 to max_W)
    W = (frame / 50) * max_W

    # Create new model with current W
    model = Square_Hamiltonian(length=length, t=t, W=W, phi=0.0, q=q)

    # Collect data for all p/q values
    new_phis = []
    new_energies = []

    for q_val in range(1, model.max_q + 1):
      for p in range(q_val + 1):
        if gcd(p, q_val) == 1:
          model.phi = p / q_val
          model.construct_hamiltonian()
          new_phis.extend([model.phi] * model.N)
          new_energies.extend(model.evals.tolist())

    # Update scatter plot
    scatter.set_offsets(np.c_[new_phis, new_energies])
    title.set_text(f'Hofstadter Butterfly (W = {W:.2f})')

    return scatter, title

  anim = FuncAnimation(
      fig,
      update,
      frames=51,  # 0 to max_W in 50 steps
      interval=100,
      blit=True
  )

  if save_path:
    print(f"Saving animation to {save_path}")
    anim.save(save_path, writer='pillow', fps=10)
  else:
    plt.show()

  plt.close()

animate_butterfly_evolution(save_path='butterfly_evolution.gif')
