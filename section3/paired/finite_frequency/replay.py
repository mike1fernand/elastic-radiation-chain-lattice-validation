"""Replay the two pointwise enclosures and their rigorously separated outputs."""
import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

if not __debug__:
    raise SystemExit("Run Python with assertions enabled (omit -O or -OO).")

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--independent-interval", action="store_true",
                    help="also run an independent 60-digit interval inverse (requires mpmath)")
args = parser.parse_args()
root = Path(__file__).resolve().parent
upstream = root / "upstream"


def read_json(file):
    return json.loads(file.read_text(encoding="utf-8"))


def run(script, *args):
    return subprocess.run([sys.executable, str(script), *map(str, args)], check=True,
                          capture_output=True, text=True).stdout


# Bind every upstream file to its recorded content and Git blob.
sources = read_json(upstream / "sources.json")
assert sources["commit"] == "6eac4ca3e89195a6b92445cfd408d9e1a8b1e26c"
for name, entry in sources["files"].items():
    data = (upstream / name).read_bytes()
    assert hashlib.sha256(data).hexdigest() == entry['sha256'], name
    header = f"blob {len(data)}\0".encode("ascii")
    assert hashlib.sha1(header + data).hexdigest() == entry['git_blob'], name

chain = read_json(upstream / "VALIDATED_PROOF_CHAIN.json")
binding = read_json(upstream / "KERNEL_CENTER_BINDING_CHECK.json")
for name, (tag, expected) in {
    "G": ("rounded_up_eG", Fraction(1,250)), "Q": ("rounded_up_eQ", Fraction(1,500))
}.items():
    exact = Fraction(binding[f'e{name}_replayed_exact_rational'])
    assert 0 < exact <= expected
    assert Fraction(binding['safe_caps'][f'e{name}']) == expected
    assert Fraction(chain['kernel_caps'][tag]) == expected
    assert Fraction(chain['kernel_caps'][f'e{name}_exact_saved_bound']) == exact
assert binding["hashes"]["kernel"] == sources["files"]["REFERENCE_KERNEL.npz"]["sha256"]
assert binding["G_exact_center_matches_targets"] and binding["safe_caps_cover_center_shift"]
assert Fraction(binding["Q_float_sum_to_certified_target_max_entry"]) == 0
assert all(binding["saved_evidence_status"].values())
assert chain["upstream_status"]["combined_replay_with_physical_and_bridge_status_certified"]
# These assertions bind existing inputs. This local replay covers the finite response,
# not the 222 upstream contour jobs.
for file, role in (
    ("KERNEL_CENTER_BINDING_CHECK.json", "exact_center_binding_and_upward_caps"),
    ("verify_point_rational.py", "point_verifier"),
    ("POINT_RATIONAL_CERTIFICATE.json", "point_receipt"),
    ("POINT_RESOLVENT_PROOF.md", "point_proof")
):
    assert sources["files"][file]["sha256"] == chain["files"][role]["sha256"]

# Run a hash-verified copy of the original ramp verifier. Redirect only
# its historical kernel path. Write the new receipt to a temporary directory.
with TemporaryDirectory(prefix="paper1_ramp_replay_") as td:
    temp = Path(td)
    source = (upstream / "verify_point_rational.py").read_text(encoding="utf-8")
    source, count = re.subn(r"^SOURCE=Path\(.*\)$", "SOURCE=HERE/'REFERENCE_KERNEL.npz'",
                             source, flags=re.MULTILINE)
    assert count == 1
    runner = temp / "verify_ramp.py"
    runner.write_text(source, encoding="utf-8")
    shutil.copyfile(upstream / "REFERENCE_KERNEL.npz", temp / "REFERENCE_KERNEL.npz")
    run(runner)
    replay = read_json(temp / "POINT_RATIONAL_CERTIFICATE.json")
pinned = read_json(upstream / "POINT_RATIONAL_CERTIFICATE.json")
assert replay["source_kernel_sha256"] == binding["hashes"]["kernel"]
for key in ("status","physical_parameters","operator_caps","tighter_conditional_interval"):
    assert replay[key] == pinned[key]
nominal = replay["tighter_conditional_interval"]
assert nominal == chain["nominal_infinite_lattice_interval"]
ramp_lo = Fraction(nominal["lower_rational"])
ramp_hi = Fraction(nominal["upper_rational"])

swap = json.loads(run(root / "verify_swapped_nominal.py", upstream / "REFERENCE_KERNEL.npz"))
expected_swap = read_json(root / "swapped_result.json")
assert swap == expected_swap
assert swap['kernel_sha256'] == replay["source_kernel_sha256"]
swap_lo = Fraction(swap["qT_lower"])
swap_hi = Fraction(swap["qT_upper"])
gap = ramp_lo - swap_hi
assert 0 < swap_lo < swap_hi < ramp_lo < ramp_hi < 1
assert gap > Fraction(49,1000)
print("PASSED: both contact responses use the same archived bare-kernel bounds.")
print(f"Ramp nominal interval         [{ramp_lo}, {ramp_hi}]")
print(f"Swapped nominal interval      [{swap_lo}, {swap_hi}]")
print(f"Difference lower bound         {gap} of incident power")
print(f"Approximate lower bound        {100*float(gap):.7f} percentage points")
print("Strict conservative result     greater than 4.9 percentage points")
print(run(root / "verify_numerical_estimates.py"),end="")
if args.independent_interval:
    print(run(root / "verify_interval_independent.py"),end="")
