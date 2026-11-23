# Low-frequency check for the ramp and the interior exchange

This calculation uses the same two six-segment contacts as the finite-frequency comparison at angular frequency 1.2. It checks both leading coefficients at the three dimensionless frequencies 0.2, 0.05 and 0.0125. No bulk integration is required to replay the check.

## Inputs and outputs

The axial line has identical leads with stiffness 1.69. In reference units the spacing, line and bulk masses, and axial and face-diagonal bulk spring stiffnesses all equal one. The sole incident wave enters from the left. From that end the six per-link stiffnesses are (1,2,3,4,5,6)/7 for the ramp and (1,2,5,4,3,6)/7 for the exchanged profile. One strength applies to all eight contact links on its segment. The reference taper (0.1,0.2,0.2,0.2,0.2,0.1) is replayed solely as an independent check against the archived data.

T_qT and T_qL are the time-averaged powers radiated into the lower two bulk branches combined and the upper branch, respectively, each divided by the incident line power. Each plotted or tabulated ratio divides one of these fractions by its own C times frequency squared. In the reference units C_T = 2(1.3)/(3π) and C_L = 1.3/[3π(√3)³]. Values near one test the leading low-frequency prediction.

## Reproduce

Extract the ZIP and run the following commands in the directory containing replay.py, using Python 3 and NumPy. All inputs are included.

    python replay.py
    python replay.py --json

The second command prints full-precision results, including the four-channel energy balance. numeric_receipt.json records the original run. Machine-dependent last digits can vary. The driver checks each included upstream file against its Git blob and SHA-256 in provenance.json, calls the retained scatter function for all profiles, and independently assembles the 48 inclined springs and solves the 19-by-19 line-and-force equations as a single system. Numerically independent quadrature resolutions use the same archived integration procedure (24 angular and 96 radial points in its size parameters, versus 40 and 128). All taper results are checked against the matching upstream rows at unit scale.

## Interpretation

| ω | Ramp qT | Ramp qL | Exchange qT | Exchange qL |
|---|---:|---:|---:|---:|
| 0.2 | 1.02720924 | 1.03115517 | 1.02395764 | 1.02820898 |
| 0.05 | 1.00167206 | 1.00192134 | 1.00147879 | 1.00174574 |
| 0.0125 | 1.00010439 | 1.00011998 | 1.00009235 | 1.00010904 |

Entries are ratios to the respective predicted limits as defined above. Both archived resolutions gave the displayed values on the verification run. These three samples illustrate convergence for these two fixed profiles. The uniform result for all admissible profiles follows from the analytic theorem. This numerical replay supplies no rigorous quadrature error bounds, onset frequency, or guarantee between the sampled frequencies.

## Provenance

The solver, kernel files, status records and reference-taper rows in upstream/ come from the [low-frequency directory at the pinned commit](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/tree/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/predictions/baseline/Prediction_Calculations_2026-09-26/low_frequency). provenance.json provides the exact path and checksums of every included upstream file. The ramp and exchanged-profile low-frequency outputs in this ZIP are new evaluations using that solver and those data. The finite-frequency enclosure at 1.2 and its separate verifier belong to the interior-swap verification package.
