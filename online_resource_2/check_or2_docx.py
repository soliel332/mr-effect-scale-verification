#!/usr/bin/env python3
"""Check the Word rendering of Online Resource 2 against the exported workbook, cell by cell.

usage: python online_resource_2/check_or2_docx.py OR2_DOCX [OR2_XLSX]

The Word file rounds for display: beta and SE to four significant figures (trailing zeros kept),
P to two significant figures in x10 notation, allele frequency to three decimals, F to one
decimal, positions with thousands separators. Sample-size columns are omitted from the Word
tables and the proxy columns appear only for FinnGen outcomes. This script applies those rules
to the workbook written by export_or2.py and requires every Word cell, header and row order to
match. It needs python-docx and openpyxl.
"""
import os
import sys

import pandas as pd
from docx import Document

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCX = sys.argv[1]
XLSX = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "online_resource_2", "Online_Resource_2_instrument_level_data.xlsx")
SUP = str.maketrans("-0123456789", "\u207b\u2070\u00b9\u00b2\u00b3\u2074\u2075\u2076\u2077\u2078\u2079")

BASE = ["SNP", "Chr", "Position", "EA", "OA", "EAF (exposure)", "EAF (outcome)", "Beta (exposure)",
        "SE (exposure)", "P (exposure)", "Beta (outcome)", "SE (outcome)", "P (outcome)", "F", "Palindromic"]
PROXY = ["Proxy", "Proxy SNP"]
INDEX = ["Table", "Exposure", "Exposure ID", "Outcome", "Outcome ID", "Variants", "Palindromic", "Outcome proxies"]


def p_fmt(v):
    s = f"{v:.2g}"
    if "e" in s:
        m, e = s.split("e")
        return f"{m}\u00d710{str(int(e)).translate(SUP)}"
    return s


def sig4(v):
    return f"{v:#.4g}"


def eaf(v):
    return v if isinstance(v, str) else f"{v:.3f}"


def rows_for(x, finngen):
    out = []
    for r in x.itertuples(index=False):
        row = [r[0], str(r[1]), f"{r[2]:,}", r[3], r[4], eaf(r[5]), eaf(r[6]), sig4(r[7]), sig4(r[8]),
               p_fmt(r[9]), sig4(r[11]), sig4(r[12]), p_fmt(r[13]), f"{r[15]:.1f}", r[16]]
        if finngen:
            row += [r[17], r[18]]
        out.append(row)
    return out


def main():
    doc = Document(DOCX)
    tables = [[[c.text for c in row.cells] for row in t.rows] for t in doc.tables]
    xl = pd.read_excel(XLSX, sheet_name=None)
    idx = xl["Index"].dropna(subset=["Sheet"])
    problems = []
    want_index = [INDEX] + [[r["Table"], r["Exposure"], r["Exposure ID"], r["Outcome"], r["Outcome ID"],
                             str(int(r["Variants"])), str(int(r["Palindromic"])), str(int(r["Outcome proxies"]))]
                            for _, r in idx.iterrows()]
    if tables[0] != want_index:
        problems.append("index table differs")
    for k, (_, r) in enumerate(idx.iterrows(), start=1):
        finngen = r["Outcome ID"] == "finn-b-J10_ASTHMA"
        want = [BASE + (PROXY if finngen else [])] + rows_for(xl[r["Sheet"]], finngen)
        got = tables[k]
        if got != want:
            diff = [(i, j, g, w) for i, (gr, wr) in enumerate(zip(got, want)) for j, (g, w) in enumerate(zip(gr, wr)) if g != w]
            problems.append(f"{r['Table']}: {len(got)} vs {len(want)} rows; {len(diff)} differing cells; first {diff[:3]}")
    n_cells = sum(len(row) for t in tables for row in t)
    if problems:
        print("\n".join(problems))
        sys.exit(1)
    print(f"Word tables match the workbook: {len(tables)} tables, {n_cells} cells, labels, order and rounding identical")


if __name__ == "__main__":
    main()
