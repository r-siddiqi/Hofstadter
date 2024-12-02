# Contributing Guidelines

Thank you for your interest in contributing to Simulating Anderson Localization and Hofstadter Butterflies on Common Euclidean Lattice Structures! This document outlines the process for contributing to this project.

# Contributing Guidelines

Thank you for your interest in contributing to our project on Simulating Anderson Localization and Hofstadter Butterflies on Common Euclidean Lattice Structures! This document outlines how to contribute effectively to this scientific computing project.

## Table of Contents
- [Getting Started](#getting-started)
- [Development Process](#development-process)
- [Submission Guidelines](#submission-guidelines)
- [Project Structure](#project-structure)
- [Testing Guidelines](#testing-guidelines)
- [Documentation Standards](#documentation-standards)

## Getting Started

Our code uses standard Python scientific libraries (NumPy, Matplotlib). 

1. Fork the repository
2. Clone your fork: `git clone https://github.com/YOUR_USERNAME/repository_name`


## Development Process

1. Create a new branch for your feature:
   ```bash
   git checkout -b feature/your-feature-name
   ```

2. Make your changes following our coding standards:
   - Use clear, descriptive variable names
   - Document functions and classes using NumPy docstring format
   - Include type hints for Python functions
   - Keep functions focused and modular
   - Follow PEP 8 style guidelines

3. Write or update tests as needed
4. Update documentation if you're introducing new features or changing existing functionality

## Submission Guidelines

1. Run and test your code thoroughly
2. Update documentation if necessary:
   - Update docstrings
   - Add example usage if appropriate
   - Update README.md if needed

3. Commit your changes:
   ```bash
   git add .
   git commit -m "Your descriptive commit message"
   ```

4. Push to your fork:
   ```bash
   git push origin feature/your-feature-name
   ```

5. Create a Pull Request (PR):
   - Use a clear, descriptive title
   - Describe your changes in detail
   - Reference any related issues

## Project Structure

```
├── lattice_types/    # Different lattice implementations
├── plots/            # Visualization tools
├── scratch_notebooks/  # Working notebooks
├── utils/            # Helper functions and utilities
├── test.ipynb        # Test notebook
└── README.md         # Project documentation
```

## Testing Guidelines

- Test your code thoroughly
- If applicable, include test cases for edge cases and boundary conditions
- Verify results against known physical expectations found in the literature
- Ensure visualizations are clear and properly labeled
- Consider testing with different lattice sizes and parameters

## Documentation Standards

### Code Documentation
- Use NumPy style docstrings
- Include type hints
- Document parameters, returns, and raises
- Provide examples in docstrings when helpful

Example:
```python
def calculate_hofstadter_butterfly(
    size: int,
    flux_range: tuple[float, float],
    num_points: int
) -> tuple[np.ndarray, np.ndarray]:
    """
    Calculate the Hofstadter butterfly spectrum for a square lattice.

    Parameters
    ----------
    size : int
        Size of the square lattice
    flux_range : tuple[float, float]
        Range of magnetic flux to sweep (min, max)
    num_points : int
        Number of flux points to calculate

    Returns
    -------
    tuple[np.ndarray, np.ndarray]
        Magnetic flux values and corresponding energy eigenvalues

    Examples
    --------
    >>> flux, energies = calculate_hofstadter_butterfly(20, (0, 1), 100)
    >>> plt.scatter(flux, energies, s=1)
    """
    # Implementation
```

### Commit Messages
- Be descriptive but concise
- Reference issues where appropriate

## Questions or Issues?

Feel free to open an issue if you have questions or run into problems. We will try our best to help!
