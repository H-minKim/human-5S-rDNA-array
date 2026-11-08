# -*- coding: utf-8 -*-
r"""Records the files a script actually read, and runs the original command unchanged (2026-10-21, for the code repository listing, figure output is identical)
   Usage: python _tracewrap.py <log.tsv> <original script.py> [args ...]
   Collects the paths of files opened in read mode from the 'open' event of sys.addaudithook (writes are marked separately)."""
import os, sys, runpy, atexit
OUTT = sys.argv[1]
SEEN = {}


def _hook(ev, args):
    if ev == "open" and args and isinstance(args[0], (str, bytes, os.PathLike)):
        p = os.fsdecode(args[0]) if not isinstance(args[0], str) else args[0]
        mode = args[1] if len(args) > 1 and isinstance(args[1], str) else "r"
        k = "w" if any(c in (mode or "") for c in "wax+") else "r"
        try: p = os.path.abspath(p)
        except Exception: pass
        SEEN.setdefault(p, set()).add(k)


def _dump():
    with open(OUTT, "w", encoding="utf8", newline="\n") as g:
        for p, ks in sorted(SEEN.items()):
            g.write("%s\t%s\n" % ("".join(sorted(ks)), p))


atexit.register(_dump)
sys.addaudithook(_hook)
sys.argv = sys.argv[2:]
sys.path.insert(0, os.path.dirname(os.path.abspath(sys.argv[0])))
runpy.run_path(sys.argv[0], run_name="__main__")
