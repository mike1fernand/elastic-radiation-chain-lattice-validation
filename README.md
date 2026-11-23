# Elastic Radiation from a Chain into a Three-Dimensional Lattice

Code and data for Mike Fernandes’s paper. A compression wave travels along a chain and sends power into a surrounding lattice through a set of springs. The paper asks whether moving stiffness between segments changes how much power goes into shear-like and longitudinal-like waves. This repository checks the contact forces, the low-frequency radiation law, and both numerical tables.

## Run the checks

1. Extract the repository ZIP or clone the repository. Open a terminal in the folder that contains `validate.py` and `requirements.txt`. Leave the companion ZIP intact. Its files are already unpacked under `section3/`.
2. Use Python 3.12 and run the commands below for your operating system. They create a local Python environment, install the two dependencies, and run the checks. Installation normally needs internet access. After that, the last command runs offline, including from the downloaded ZIP.

**Windows (PowerShell)**

```powershell
py -3.12 -m venv .venv
.\.venv\Scripts\python.exe -m pip install -r requirements.txt
.\.venv\Scripts\python.exe validate.py
```

**macOS or Linux**

```sh
python3.12 -m venv .venv
.venv/bin/python -m pip install -r requirements.txt
.venv/bin/python validate.py
```

**How to tell it worked.** The first six tests print `OK`. Let the run continue. Near the end, look for `Passed both table comparisons, the baseline and the two nominal-response enclosures.` The last line reports the lower bound on the decrease in shear-like radiation. You can run the checks without installing LaTeX or downloading the long paper’s repository.

## What the results mean

| Paper location | What is checked |
| --- | --- |
| Section 1 and Supplement S1 | Six exact tests check the spring geometry, the forces on the chain and lattice, and the stiffness exchange. |
| Section 2 and Supplement S2 | The checks cover the lattice’s long-wavelength response and the two leading radiation coefficients. The two leading coefficients stay the same when stiffness is redistributed among fixed attachment positions at fixed positive total stiffness. Eight arrangements at three frequencies illustrate that limit. The proof for all allowed arrangements is in S2. |
| Section 3, Table 1 | For the two arrangements, the calculation reproduces the percentage by which each radiated fraction differs from its low-frequency prediction. |
| Section 3, Table 2 | With the wave entering from the same end at dimensionless frequency 1.2, exchanging two stiffnesses lowers the estimated shear-like radiation from **73.6% to 67.9% of incident chain power**. The estimated decrease is **5.7 percentage points**. Using the archived error bounds on the lattice response, the decrease exceeds **4.9 percentage points** for these two contacts at this frequency. |

## Where to go next

- [Check the tables against your manuscript or save the unrounded results](paper/README.md).
- [See how the forces and low-frequency coefficients are checked](supplement_checks/README.md).
- [Find the original lattice calculations in the long paper’s repository](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/tree/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c). This repository includes the inputs used in the shorter paper.
