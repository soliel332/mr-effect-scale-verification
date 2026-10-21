#!/usr/bin/env python3
"""Run the code example in README.md ('Applying the check to another release') and check its printed result."""
import contextlib, io, os, re, sys

ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
readme = open(os.path.join(ROOT, "README.md"), encoding="utf-8").read()
section = readme.split("## Applying the check to another release", 1)[1].split("\n## ", 1)[0]
code = re.search(r"```python\n(.*?)```", section, re.S).group(1)
expected = re.search(r"^# (.+)$", code, re.M).group(1).strip()
os.chdir(ROOT)
buf = io.StringIO()
with contextlib.redirect_stdout(buf):
    exec(compile(code, "README example", "exec"), {})
got = buf.getvalue().strip()
assert got == expected, (got, expected)
print(f"README example reproduced: {got}")
