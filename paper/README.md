# Reproducing Tables 1 and 2

The main validation command checks [the two tables extracted from the manuscript](tables.tex). It recalculates their values and checks the column headings and the parts of each caption that define what the numbers mean. To check a revised manuscript, give the program your edited `Manuscript.tex` from the separate paper build.

## Check an edited manuscript

Start in the repository folder containing `validate.py` after completing [the Python setup](../README.md#run-the-checks). Replace the example path with your file’s path. Keep the quotation marks if a folder name contains spaces.

**Windows (PowerShell)**

```powershell
.\.venv\Scripts\python.exe validate.py --tex "C:\path to paper\Manuscript.tex"
```

**macOS or Linux**

```sh
.venv/bin/python validate.py --tex "/path to paper/Manuscript.tex"
```

Keep the LaTeX comments that mark the start and end of each table’s caption, header and rows. For example, the first caption sits between `% TABLE1-CAPTION-BEGIN` and `% TABLE1-CAPTION-END`. The corresponding `HEADER` and `ROWS` pairs are used too, as are all three pairs for Table 2.

## See or save the unrounded results

The [original numerical record](../section3/receipt.json) includes the four power fractions, the two tables and the bounds for the two contacts. To write a fresh record, run the command for your system from the repository folder. It creates or replaces `local_receipt.json` there.

**Windows (PowerShell)**

```powershell
.\.venv\Scripts\python.exe section3/check_tables.py --tex paper/tables.tex --output local_receipt.json
```

**macOS or Linux**

```sh
.venv/bin/python section3/check_tables.py --tex paper/tables.tex --output local_receipt.json
```

In Table 1, a value of 3.1 means that the computed fraction is 3.1% above its own low-frequency prediction. In Table 2, each value is a percentage of the **incident chain power**. Table 2 uses surface integration. For the exchanged contact, that method and tetrahedron integration both round to 67.9% of incident power. The two independent grid calculations round to 67.8%. The decrease rounds to 5.7 percentage points with all four methods when calculated from each method’s unrounded results. The lower bound of 4.9 percentage points uses the archived errors on the lattice response and applies to these two arrangements at the tested frequency.

## Follow the calculation back to its sources

Supplement S3 cites [the original companion ZIP](../Section_3_Table_Validation.zip). All 41 files are already in [section3/](../section3/), and the validation program checks them against that ZIP. Its SHA-256 fingerprint is `32ddbffabaf6e82637ee82c17190c54137b5b57627aaca5a322da1b8ecd7818e`. The manuscript from which [tables.tex](tables.tex) was extracted has SHA-256 fingerprint `68c673e59f514c0f030f02ccc8e1197f4fa01d6dda9ee968737f6af448bb8518`.

The original [companion guide](../section3/README.md) and [guide to the two contacts](../section3/paired/README.md) retain the instructions shipped with the archive. Their “outline” was an earlier draft. Use the manuscript command above for the current paper. In the shorter paper, Section 2 and Supplement S2 give the low-frequency theorem and proof. Section 3 and Supplement S3 give the numerical comparison. The archived [low-frequency calculation](../section3/paired/low_frequency/README.md) and [finite-frequency bound](../section3/paired/finite_frequency/README.md) describe how the saved inputs were used.
