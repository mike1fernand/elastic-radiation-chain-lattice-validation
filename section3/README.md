# Section 3 table validation

The two tables concern the same six-segment contact with the wave incident from the left. The actual stiffness of each of the eight springs at segment j is its listed number times K/7. The ramp is (1, 2, 3, 4, 5, 6) and the exchange is (1, 2, 5, 4, 3, 6). Both use equal bulk spring stiffnesses K, identical line leads, equal site masses and chain spring stiffness 1.69 K. Angular frequency Ω is in units of √(K/m_b).

## Replay

Extract the directory and run from `section3/` with Python 3.10 or later, NumPy and mpmath. The supplied receipt used CPython 3.12.14, NumPy 2.3.5 and mpmath 1.3.0. Keep assertions enabled. The calculation reads the included files and needs no network access.

```sh
python check_tables.py
python check_tables.py --tex /path/to/section3_outline.tex --output local_receipt.json
```

The first command checks the marked caption, header and body of each table in `table_rows.tex`. The caption checks cover the signed percentage-departure formula for Table 1 and the power denominator, percentage units, test frequency and estimate method for Table 2. Keep the `TABLE1-CAPTION` and `TABLE2-CAPTION` markers around the visible captions when transferring the tables into the manuscript. The second checks the identically marked tables in the supplied LaTeX source and writes all unrounded channel fractions, low-frequency diagnostics, derived fractions, displayed strings and exact nominal endpoints. `receipt.json` contains these results for the supplied rows. The last digits of floating-point results can vary with hardware and numerical libraries. The rational endpoints are exact. As a negative control, copy the table source, replace Table 2’s “percentages of incident chain power” with “percentages of emitted bulk power”, and pass the modified copy to `--tex`. It must stop with a power-denominator error.

## What is checked

- At Ω = 0.2, 0.05 and 0.0125, the archived harmonic solver and a second assembly from the 48 springs check four output fractions at two quadrature resolutions. At each resolution the code also compares the original reference taper against its archived unit-strength rows as a control. Table 1 displays the signed departure `100 × [T_B/(C_B ω²) − 1]` to two significant digits for each branch group. The program converts each group’s **own** coefficient to the stated reference units and checks the displayed cells at both resolutions. The three samples illustrate the theorem of Section 2.
- Table 2 uses the archived surface integration and percentages of incident chain power, in the order reflected, onward, combined lower pair and upper branch. Both contacts are solved with the same left incidence. The program checks the four unrounded entries against their pinned record, verifies four-channel closure, and separately recomputes the tetrahedron and two independent grid estimates. The exchange shear-like estimate is 67.8571044% from the surface kernel, 67.8501620% from the tetrahedron kernel and 67.8499104% from the independent N = 768 grid. Hence the first two round to 67.9%. The grid rounds to 67.8%. Subtract the two unrounded incident fractions before converting their difference to percentage points. The receipt also identifies the share with only emitted bulk power in its denominator.
- The paired finite replay checks the original rational ramp receipt, the separate exchange certificate, their common bare-lattice data and the 60-digit interval implementation. Subtracting the exchange upper endpoint from the ramp lower endpoint gives a lower bound of 4.9853073 percentage points. Report a decrease **greater than 4.9 percentage points** for these two nominal contacts at Ω = 1.2. The stored upstream contour certificate and its target binding supply the operator error bounds on the uncoupled lattice. The small replay propagates them through the two contact solutions. The numerical estimates and power balance do not certify the printed tenths or a frequency band.

## Source records

The two archived source packages and the hashes of the included source files are identified in `inputs.json`. `paired/` is the unchanged extracted `Paper1_Same_Contact_Pair_Verification.zip`. It provides both same-pair low-frequency rows and the exchange enclosure. Its finite source inventory maps the original ramp certificate and all four original numerical channels for **both** profiles to Git commit `6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c`. The center-binding receipt and proof chain identify the source of the scaled operator error bounds, 1/250 and 1/500. The full upstream source is [available at the pinned commit](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/tree/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c). The direct component entry points can also be run individually.

```sh
python paired/low_frequency/replay.py --json
python paired/finite_frequency/replay.py --independent-interval
```

`independent/` contains the two finite-frequency grids and their producer and contact solver from `Paper1_Focused_Independent_Check_2026-11-20.zip`. To regenerate just these comparison grids into separate files, use the reconstruction settings below. The independent source archive documented the N = 768 settings. For N = 512, 300 polar points reproduce the saved surface matrix to floating-point roundoff. Its historical polar count is unrecorded. Neither regeneration is needed to replay the tables.

```sh
python independent/scripts/make_green.py 1.2 512 0.03,0.05,0.07,0.09,0.11,0.13 300 new_N512.npz
python independent/scripts/make_green.py 1.2 768 0.02,0.04,0.06,0.08,0.10,0.12,0.14 400 new_N768.npz
```
