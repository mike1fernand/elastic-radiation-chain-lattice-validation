"""Reproduce the calculations for Elastic Radiation from a Chain into a Three-Dimensional Lattice."""
import argparse
import hashlib
from pathlib import Path
import subprocess
import sys
from zipfile import ZipFile

ROOT = Path(__file__).resolve().parent
COMPANION_HASH = "32ddbffabaf6e82637ee82c17190c54137b5b57627aaca5a322da1b8ecd7818e"


def verify_companion():
    archive = ROOT / "Section_3_Table_Validation.zip"
    if hashlib.sha256(archive.read_bytes()).hexdigest() != COMPANION_HASH:
        raise ValueError("The companion ZIP differs from the archive cited in Supplement S3.")
    with ZipFile(archive) as bundle:
        for member in bundle.infolist():
            if member.is_dir():
                continue
            if (ROOT / member.filename).read_bytes() != bundle.read(member):
                raise ValueError(f"{member.filename} differs from the cited ZIP.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--tex", type=Path, default=ROOT / "paper" / "tables.tex",
                        help="Check the two marked tables in this manuscript or excerpt.")
    args = parser.parse_args()
    if not __debug__:
        parser.error("Use Python with assertions enabled.")
    tex = args.tex.resolve(strict=True)
    verify_companion()
    print("Confirmed the companion archive and its extracted files.", flush=True)
    for script in ["check_section1.py", "check_section2.py"]:
        subprocess.run([sys.executable, str(ROOT / "supplement_checks" / script)],
                       cwd=ROOT, check=True)
    subprocess.run([sys.executable, str(ROOT / "section3" / "check_tables.py"), "--tex", str(tex)],
                   cwd=ROOT, check=True)


if __name__ == "__main__":
    try:
        main()
    except (ValueError, OSError, subprocess.CalledProcessError) as error:
        sys.exit(str(error))
