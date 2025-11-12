---
file_format: mystnb
kernelspec:
  name: python3
---

```{image} _static/logo.png
:align: center
:alt: ScopeSim-Targets logo
:width: 400 px
```

# ScopeSim Targets
Simplified, unified target description model for the `ScopeSim` observation simulation engine.

## Getting started
### Installation
Recommended way to install the package is to simply use `pip`:

```bash
pip install scopesim-targets
```

### Basic usage
In the most basic example, we will set up a simple point source and observe it with `ScopeSim` using METIS imaging mode:

```{code-cell} ipython3
:tags: [hide-output]

from scopesim import Simulation
from scopesim_targets.point_source import Star

target = Star(
    position=(0, 0),
    spectrum="G2V",
    brightness=("K", 18)
)

simulation = Simulation("MICADO")
simulation(target.to_source())
```

```{code-cell} ipython3
img_slice = slice(462, 562)  # only show center of image
simulation.plot(img_slice=(img_slice, img_slice))
```

## Contents

```{toctree}
:maxdepth: 2

defining_positions.md
defining_spectra.md
defining_brightness.md
yaml_syntax.md
```

## API reference

```{eval-rst}
.. autosummary::
   :toctree: _autosummary
   :template: custom-module-template.rst
   :recursive:
   :caption: Package Contents

   scopesim_targets.target
   scopesim_targets.point_source
```
