# Verification for the same two contacts

The two folders test the same original and exchanged contacts, from the same incident side. The third and fifth of six per-segment stiffnesses are exchanged. Run the examples after extracting this entire ZIP. The manuscript outline gives the intended role of each result.

## Low-frequency numerical check

In the low_frequency folder run python replay.py (Python 3 and NumPy). The program checks the two bulk power fractions against their separate leading laws at three frequencies. Each recorded power is divided by incident line power, then by its own coefficient times frequency squared in the displayed ratios. This relies on archived quadratures. The uniform limit is an analytical result in Paper 1, Section 5.3 and Supplement SM7A.

## Finite-frequency change

In the finite_frequency folder run python replay.py (Python 3 and NumPy). To repeat the independent interval calculation, also install mpmath and run python replay.py --independent-interval from that folder. This checks both response enclosures at dimensionless angular frequency 1.2 using the error bounds already established for the bulk matrices. The difference between the estimated shear-like fractions is about 5.73 percentage points. Under the archived bounds the drop exceeds 4.9 percentage points.

Each folder contains its own README, code, data, verification records and links to the upstream files at commit 6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c. The same-pair low-frequency rows and the finite-frequency enclosure for the exchange were calculated from those inputs. The expensive bulk integrals were not recomputed for this outline.
