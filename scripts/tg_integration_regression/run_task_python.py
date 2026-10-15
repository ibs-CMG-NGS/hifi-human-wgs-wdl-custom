#!/usr/bin/env python3
"""Run the embedded python3 heredoc of a transgene_integration task outside miniwdl.

Usage: run_task_python.py <task_name> <workdir> key=value [key=value ...]
~{key} placeholders in the task's python block are replaced by the given values.
Used for quick regression checks of extract_integration_coords / detect_chimeric_reads on saved PAFs.
"""
import os
import re
import subprocess
import sys

WDL = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "workflows", "wdl-common", "wdl", "tasks", "transgene_integration.wdl")


def task_python(task):
    s = open(WDL).read()
    i = s.index(f"task {task} {{")
    j = s.find("\ntask ", i + 1)
    block = s[i:j if j > 0 else len(s)]
    m = re.search(r"python3 << 'PYEOF'\n(.*?)\nPYEOF", block, re.S)
    return m.group(1)


if __name__ == "__main__":
    task, workdir, *kv = sys.argv[1:]
    code = task_python(task)
    for k, v in (x.split("=", 1) for x in kv):
        code = code.replace("~{" + k + "}", v)
    left = re.findall(r"~\{[^}]+\}", code)
    if left:
        sys.exit(f"unsubstituted placeholders: {sorted(set(left))}")
    os.makedirs(workdir, exist_ok=True)
    sys.exit(subprocess.run([sys.executable, "-"], input=code, text=True, cwd=workdir).returncode)
