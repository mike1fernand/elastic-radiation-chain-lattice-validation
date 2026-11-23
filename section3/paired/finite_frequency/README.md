# Finite-frequency check for the exchange of two interior segment stiffnesses

This calculation bounds the fraction of incident line power emitted into the two lower bulk branches at dimensionless frequency 6/5. The six motif strengths from the incident end are (1,2,5,4,3,6)/7. The reference ramp is (1,2,3,4,5,6)/7. The check uses the same six positions, bare bulk, reference units, line stiffness 169/100 and left incident wave as the archived ramp calculation. Each strength is shared by the eight springs of its motif. At this frequency the swap reduces the exact shear-like fraction by more than 4.9 percentage points, using the bare-lattice kernel enclosures supplied in the pinned validation repository.

## Reproduce the numerical enclosure

Extract the ZIP into any writable directory. Open a terminal in the extracted directory containing replay.py. Python 3 and NumPy are required. The inputs are included. Run with assertions enabled.

    python replay.py

The script checks the included upstream files against their recorded Git blobs and SHA-256 hashes. It replays the original nominal ramp verifier with the kernel path redirected to the included file, then runs the new verifier for the swapped contact. It compares exact rational output endpoints, verifies their strict separation, and recomputes the numerical estimates using both recorded kernels. It creates the original ramp verifier’s receipt in a temporary directory.

To also check the two contact systems using a separately implemented 60-decimal-digit interval-matrix inversion, install mpmath and run:

    python replay.py --independent-interval

The new swapped-contact verifier uses a floating inverse as a proposed preconditioner. Subsequent Fraction comparisons certify its residual, metric positivity, operator bounds and final enclosure. Rational brackets for the two irrational model constants enclose sqrt(2) and sqrt(133). sources.json in the upstream folder records each upstream file and its path in the pinned repository. The local additions include the replay driver, the swapped-contact verifier and receipt, and two cross-checks.

## Formula and proof dependencies

The pinned source is [the nominal ramp verifier](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/blob/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/ramp/tolerance/verify_point_rational.py), SHA-256 8cc7821ddc6c2adb20f9ff3079cda17519086479e11b1ddeb833053dcdcf3359. Its [proof](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/blob/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/ramp/tolerance/POINT_RESOLVENT_PROOF.md) gives the residual and source-to-output inequalities. This verifier changes the third and fifth strengths and checks new bounds on all quantities that depend on the contact. It takes the kernel path as an argument.

Both cases use the same 12-dimensional bare-lattice kernel center SHA-256 82777dacc866367cc2a4ffcb9b5430698afdcf577984a10186a48ab61f28c52b. The operator error caps eg = 1/250 for the scaled retarded Green matrix and eq = 1/500 for its scaled shear-like spectral matrix come from [the existing center-binding record](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/blob/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/ramp/bound_review/KERNEL_CENTER_BINDING_CHECK.json). They concern the bare medium at fixed frequency on the six fixed positions. Exchanging two segment stiffnesses changes the finite response and leaves these bare input operators the same. The full-domain contour certificate is inherited. This small verifier does not recompute it.

For the swapped profile let h bound the source resolvent norm, j the response to output norm, and z0 the nominal source norm, in the normalization of the pinned proof. The script checks h = 0.95, j = 0.642, z0 = 0.57, and nominal output amplitude in [0.82375, 0.82376]. Set Z = z0/(1-h eg) and d = j eg Z. The exact fraction lies between (0.82375-d)^2 - 2 eq Z^2 and (0.82376+d)^2 + 2 eq Z^2. The exact dyadic metric has a verified rational LDL factorization with positive pivots. The exact spectral metric inherits positivity from its spectral definition.

The resulting swapped-contact endpoints are 6697169359528801/9924144400000000 and 661266871782241/969154726562500. They lie strictly between 0.674 and 0.683. The pinned [nominal ramp proof chain](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/blob/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/ramp/tolerance/VALIDATED_PROOF_CHAIN.json) gives a ramp fraction greater than 0.732. Thus exchanging the two strengths lowers the exact fraction by more than 0.049, or 4.9 percentage points, for this single frequency and these two profiles. This does not certify the displayed tenths of a percent or a frequency interval. The [archived finite estimates](https://github.com/mike1fernand/low-frequency-radiation-reciprocal-junction-validation/blob/6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c/validation/referee_checks/profile_reversal.json) give about 73.6% and 67.9%.
