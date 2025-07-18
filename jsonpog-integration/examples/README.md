# Examples

Basic examples for standard usage of JSON tables with `correctionlib` are provided.

All examples are provided at least in Python, but some are also provided in C++. Both essentially do the same, with minor differences. The Python example can be directly executed in the command line, whereas the C++ example should first be compiled with `make`.

## BTV

Four use cases are shown
- fixedWP correction with mujets (here medium WP)
- fixedWP correction uncertainty (here tight WP and comb SF)
- shape correction SF 
- shape correction SF uncertainties

The results from C++ and Python examples differ in values, as the Python (C++) executable relies on NumPy / ROOT for the random number generation.

## Electrons

## JERC

Six use cases are shown:
- Retrieve the jet energy correction at a single level (e.g. `L2Relative`) for AK4 and AK8.
- Retrieve the jet energy compound correction (all levels---this is probably what you need in most cases).
- Retrieve the jet energy uncertainty (e.g. `Total`---probably also what most analyses will need).
- Retrieve the jet energy resolution scale factor
- Retrieve the jet energy resolution 
- Retrieve the jet energy correction factor from the jet energy resolution smearing

## JMAR

## MET $\phi$

## Muons

Retrieving various SFs for Run 2 UL and Run 3 campaigns.

## Photons

Retrieving various SFs for Run 2 UL and Run 3 campaigns.

## Tauons

