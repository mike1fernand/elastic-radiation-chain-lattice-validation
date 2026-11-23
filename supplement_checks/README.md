# Check the contact forces and low-frequency radiation

These checks accompany Sections 1 and 2 of *Elastic Radiation from a Chain into a Three-Dimensional Lattice*. The [main instructions](../README.md#run-the-checks) run them automatically.

## Run either check

Open a terminal in the folder containing `validate.py`. Use the environment from the main instructions. In the examples below, replace `python` with `.venv/bin/python` on macOS or Linux, or `.\.venv\Scripts\python.exe` on Windows. The first check uses Python alone. The second also uses NumPy and the supplied data.

```sh
python supplement_checks/check_section1.py
python supplement_checks/check_section2.py
```

For Section 1, look for `Ran 6 tests` and `OK`. Section 2 prints the relative errors for each wave group and reports `Endpoint and spring force sums agree in all 24 cases`. A failed check exits with an error.

## What the scripts check

**Section 1 and Supplement S1.** Two chain masses connect to four lattice corners through eight springs. The script checks the bond lengths, the force on each mass from the spring energy, and the balance of axial forces. It checks that translating all six masses together produces zero force. With the corners initially at equilibrium and the chain segment’s midpoint fixed, separating its two masses pulls the corners inward. Bringing the masses closer pushes them outward. A final check verifies that only segments 3 and 5 exchange stiffnesses. All six tests use exact fractions. [Section 1 script](check_section1.py)

**Section 2 and Supplement S2.** The script checks the long-wavelength lattice equation (S11) from the nine spring directions using exact fractions. It also checks the angular weights and leading radiation coefficients for the shear-like (`qT`) and longitudinal-like (`qL`) groups.

At three dimensionless frequencies (0.2, 0.05 and 0.0125), it tests the ramp, its two-segment exchange and six cases with all the stiffness on one segment. For these 24 cases it checks the force balance and compares the radiated fractions and lattice responses with their low-frequency limits. Supplement S2 gives the proof for all allowed arrangements. [Section 2 script](check_section2.py)

## Input files and optional regeneration

The scripts are adapted from the earlier focused build, `Elastic_Radiation_Sections_1_2_Build.zip` (`validation/`). The low-frequency check uses [three Cartesian response files](data/), [three reduced response files](../section3/paired/low_frequency/upstream/refined/), and the [junction solver](../section3/independent/scripts/junction.py). The Cartesian files and solver match those in `Paper1_Focused_Independent_Check_2026-11-20.zip` (SHA-256 `2705df1ae2bedb9f4855ec27ecc00210dd40a93d5d3c28c4a23be2955e85bcf6`).

The routine checks use the saved matrices. To recompute one Cartesian response on a 512-by-512-by-512 grid, run this command from the repository root. Choose a new output filename to keep earlier results.

```sh
python section3/independent/scripts/make_green.py 0.2 512 "0.03,0.05,0.07,0.09,0.11,0.13" 240 regenerated_0.2.npz
```

The arguments are the dimensionless frequency, grid width, six positive regulator values, number of polar angles and output filename. For 0.05 or 0.0125, change the frequency and filename. The earlier archive documents 240 polar angles at 0.2. For the two lower frequencies, it reports that using 240 reproduces the saved surface matrices to floating-point roundoff. The original angle counts for those two runs are unknown.
