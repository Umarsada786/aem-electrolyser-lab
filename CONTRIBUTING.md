# Contributing

Open an issue describing the scientific question, units, operating ranges and expected behaviour. For new data, include measurement provenance, permission/license, uncertainty, equipment geometry and whether records are raw measurements or derived values. Never combine simulations with experimental validation targets.

Keep calibration and evaluation sets separate. Explain any change to the split or exclusion rules. Add conservation-law or regression tests when changing model equations, and keep Python and browser implementations numerically consistent. Run `python validate.py` and `python -m unittest discover -s tests -v` before opening a pull request.

Please report fit failures, nonphysical extrapolations and browser bugs with a minimal input example. Avoid presenting good curve agreement as unique identification of material properties.
