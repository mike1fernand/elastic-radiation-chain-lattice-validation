"""Verify both Section 3 tables using the paired archive and the separate grid check."""
import argparse
import contextlib
from fractions import Fraction
from hashlib import sha256
import io
import json
from pathlib import Path
import re
import runpy
import subprocess
import sys

import numpy as np


ROOT = Path(__file__).resolve().parent
PAIRED = ROOT / "paired"
FINITE = PAIRED / "finite_frequency"
POWER_FIELDS = ("R", "Tline", "TqT", "TqL")
CONTACTS = ("ramp", "exchange")
FREQUENCIES = (0.2, 0.05, 0.0125)
LOW_COLUMNS = (("ramp", "qT"), ("swap", "qT"), ("ramp", "qL"), ("swap", "qL"))
LOW_HEADER = [r"$\Omega$", "Ramp qT", "Exchange qT", "Ramp qL", "Exchange qL"]
HIGH_HEADER = ["Contact", "Reflected", "Onward", "Shear-like pair", "Longitudinal-like"]


def require(condition, message):
    if not condition:
        raise ValueError(message)


def load(path):
    return json.loads(path.read_text(encoding="utf-8"))


def verify_inputs():
    sources = load(ROOT / "inputs.json")
    for name, digest in sources["sha256"].items():
        require(sha256((ROOT / name).read_bytes()).hexdigest() == digest, f"Changed input: {name}")


def marked_lines(tex, label):
    lines = [line.strip() for line in tex.splitlines()]
    start, end = f"% {label}-BEGIN", f"% {label}-END"
    require(lines.count(start) == lines.count(end) == 1, f"Missing or repeated markers for {label}")
    first, last = lines.index(start), lines.index(end)
    require(first < last, f"Reversed markers for {label}")
    return lines[first + 1:last]


def check_caption(tex, number):
    # Work only from the visible caption between the two markers.
    # Preserve escaped percentages and discard LaTeX line comments.
    lines = marked_lines(tex, f"TABLE{number}-CAPTION")
    caption = " ".join(re.split(r"(?<!\\)%", line, maxsplit=1)[0] for line in lines)
    prose = " ".join(caption.casefold().replace("~", " ").split())
    formulas = re.findall(r"\$([^$]+)\$|\\\((.*?)\\\)", caption)
    math = {"".join((dollar or parens).split()) for dollar, parens in formulas}
    percent = re.search(r"\bpercent(?:ages?)?\b", prose) or r"\%" in prose
    require(percent, f"Table {number} caption must identify percentage units")
    if number == 1:
        require(r"100[T_B/(C_B\omega^2)-1]" in math,
                "Table 1 caption must give the signed departure from each group's own leading law")
    else:
        require(r"\Omega=1.2" in math, "Table 2 caption must identify the test frequency")
        denominator = re.search(
            r"(?:percent(?:ages?)? of|normalized to|relative to|divided by) (?:the )?incident chain power\b",
            prose,
        )
        require(denominator, "Table 2 caption must identify the power denominator")
        require("surface integration" in prose, "Table 2 caption must identify the displayed method")
        require(re.search(r"\bestimat(?:ed|es?|ions?)\b", prose),
                "Table 2 caption must identify the entries as estimates")


def marked_rows(tex, label):
    rows = []
    for line in marked_lines(tex, label):
        require(line.endswith(r"\\"), f"Expected one LaTeX row per line in {label}: {line}")
        rows.append([value.strip() for value in line[:-2].split("&")])
    return rows


def read_table(tex, number, columns, expected):
    header = marked_rows(tex, f"TABLE{number}-HEADER")
    require(header == [columns], f"Table {number} column order or labels changed: {header}")
    actual = marked_rows(tex, f"TABLE{number}-ROWS")
    require(actual == expected, f"Table {number} mismatch\nExpected: {expected}\nFound: {actual}")


def low_table():
    cmd = [sys.executable, str(PAIRED / "low_frequency/replay.py"), "--json"]
    result = json.loads(subprocess.check_output(cmd, text=True))
    rows = {(r["kernel"], r["omega"], r["profile"]): r for r in result["rows"]}
    both = {}
    for level in ("refined", "validation"):
        display = []
        for omega in FREQUENCIES:
            line = [str(omega)]
            for profile, branch in LOW_COLUMNS:
                row = rows[level, omega, profile]
                coefficient = result["CT" if branch == "qT" else "CL"]
                ratio = row["T" + branch] / (coefficient * omega**2)
                require(np.isclose(ratio, row["ratio_" + branch], rtol=1e-14),
                        "Low-frequency power or coefficient mismatch")
                # Two significant digits, including trailing zeros.
                departure = 100 * (ratio - 1)
                cell = format(departure, "#.2g").rstrip(".")
                # The saved discrepancy is the maximum over all four channel fractions.
                uncertainty = 100 * row["second_solve_max_abs_error"] / (coefficient * omega**2)
                endpoints = (departure - uncertainty, departure + uncertainty)
                require(all(format(value, "#.2g").rstrip(".") == cell for value in endpoints),
                        f"Two contact solves can round differently at {omega} for {profile} {branch}")
                line.append(cell)
            display.append(line)
        both[level] = display
    require(both["refined"] == both["validation"], "Table 1 differs between resolutions")
    return both["refined"], result


def finite_replay():
    # Verifies upstream ramp receipt, authenticated caps, both rational intervals and strict gap.
    # Also checks the 60-digit independent interval solution against both sets of norm caps.
    cmd = [sys.executable, str(FINITE / "replay.py"), "--independent-interval"]
    subprocess.run(cmd, check=True, stdout=subprocess.PIPE, text=True)
    upstream = FINITE / "upstream"
    ref = load(upstream / "POINT_RATIONAL_CERTIFICATE.json")["tighter_conditional_interval"]
    swap = load(FINITE / "swapped_result.json")
    low, high = (Fraction(ref[key]) for key in ("lower_rational", "upper_rational"))
    swap_low, swap_high = (Fraction(swap[key]) for key in ("qT_lower", "qT_upper"))
    gap = low - swap_high
    require(0 < swap_low < swap_high < low < high < 1 and gap > Fraction(49, 1000),
            "Nominal bound or signed subtraction failed")
    return {"ramp_interval": [str(low), str(high)], "exchange_interval": [str(swap_low), str(swap_high)],
            "gap_fraction_lower_bound": str(gap), "gap_lower_bound_points": float(100 * gap)}


def checked_powers(values):
    row = {name:float(value) for name, value in zip(POWER_FIELDS, values, strict=True)}
    require(min(row.values()) >= 0 and abs(sum(row.values()) - 1) < 1e-12,
            f"Invalid four-channel powers: {row}")
    return row


def finite_tables():
    # The existing verifier also checks all four surface channels against the pinned receipt.
    with contextlib.redirect_stdout(io.StringIO()):
        calculator = runpy.run_path(str(FINITE / "verify_numerical_estimates.py"))
    ramp = np.arange(1, 7, dtype=float) / 7
    profiles = {"ramp":ramp, "exchange":ramp[[0, 1, 4, 3, 2, 5]]}
    powers = {}
    for label, key in [("surface", "surface"), ("tetrahedron", "tetra")]:
        powers[label] = {name:checked_powers(calculator["solve"](p, calculator[key], 0))
                         for name, p in profiles.items()}
    # Independent 72-coordinate calculation uses two separately stored grid matrices.
    independent_solve = runpy.run_path(str(ROOT / "independent/scripts/junction.py"))["solve"]
    for grid in (512, 768):
        with np.load(ROOT / f"independent/data/green_w1.2_N{grid}.npz", allow_pickle=False) as kernel:
            require(float(kernel["omega"]) == 1.2 and int(kernel["N"]) == grid, "Wrong independent grid")
            powers[f"grid{grid}"] = {}
            for name, profile in profiles.items():
                state = independent_solve(profile, 1.2, kernel["ReG"], kernel["AT"], kernel["AL"])
                powers[f"grid{grid}"][name] = checked_powers(state[k] for k in POWER_FIELDS)
    rounding = {label:format(100 * rows["exchange"]["TqT"], ".1f") for label, rows in powers.items()}
    require(rounding == {"surface":"67.9", "tetrahedron":"67.9", "grid512":"67.8", "grid768":"67.8"},
            f"Changed method-specific exchange rounding: {rounding}")
    display = [[name.capitalize()] + [f"{100 * powers['surface'][name][k]:.1f}\\%" for k in POWER_FIELDS]
               for name in CONTACTS]
    derived = {}
    for label, rows in powers.items():
        quantities = {}
        for name, row in rows.items():
            bulk = row["TqT"] + row["TqL"]
            quantities[name] = {"bulk_fraction_incident": bulk, "qT_fraction_of_bulk": row["TqT"]/bulk}
        quantities["ramp_minus_exchange_qT_points"] = 100*(rows["ramp"]["TqT"]-rows["exchange"]["TqT"])
        derived[label] = quantities
    return display, powers, derived


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, default=ROOT/"table_rows.tex", help="LaTeX file containing the marked table rows")
    parser.add_argument("--output", type=Path, help="Write full replay details to this JSON file")
    args = parser.parse_args()
    require(__debug__, "Run with Python assertions enabled")
    verify_inputs()
    tex = args.tex.read_text(encoding="utf-8")
    for number in (1, 2):
        check_caption(tex, number)
    table1, low = low_table()
    read_table(tex, 1, LOW_HEADER, table1)
    bounds = finite_replay()
    table2, high, derived = finite_tables()
    read_table(tex, 2, HIGH_HEADER, table2)
    if args.output:
        receipt = {"table1_cells":table1, "table2_cells":table2, "low_frequency_replay":low,
                   "finite_frequency_fractions_of_incident":high, "derived":derived, "enclosures":bounds}
        args.output.write_text(json.dumps(receipt, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print("Passed both table comparisons, the baseline and the two nominal-response enclosures.")
    for method, rows in high.items():
        print(f"{method:12} ramp {100*rows['ramp']['TqT']:.7f}%  exchange {100*rows['exchange']['TqT']:.7f}%"
              f"  change {derived[method]['ramp_minus_exchange_qT_points']:.7f} percentage points")
    print(f"Exact nominal difference exceeds 4.9 percentage points (lower endpoint {bounds['gap_lower_bound_points']:.7f}).")


if __name__ == "__main__":
    main()
