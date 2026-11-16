# -*- coding: utf-8 -*-
# A table of the basecaller release per bonobo read: a BAM MM head of "C+h?..." means guppy 6.3.7, and "C+m?..." alone means guppy 6.3.8.
#   (The same rule as x17; the v2 figures in x19 to x21 use it to exclude the 6.3.7 reads.)
#   Output: tables/BONOBO_READ_BASECALLER.tsv, with the read id, the basecaller and the MM head.
import io, struct, collections

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'Environment variables')" % k)
    return v


OUT = _need("RDNA_PREV_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
src = io.open(OUT + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
assert L11 in src
G = {"__name__": "x11hdr"}
exec(src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()"), G)
rows = {}
for refs, r in G["bam_records"](B + "/04_bam/bonobo_mPanPan1/bonobo_mPanPan1.5S.bam"):
    refid, pos0, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack("<iiBBHHHi", r[:20])
    if flag & 0x900:
        continue
    rid = r[32:32 + lrn - 1].decode()
    p = 32 + lrn + 4 * ncig + (lseq + 1) // 2 + lseq
    T = G["parse_tags"](r, p)
    if "MM" not in T:
        continue
    head = T["MM"].split(",")[0]
    rows[rid] = ("guppy-6.3.7" if T["MM"].startswith("C+h") else "guppy-6.3.8" if T["MM"].startswith("C+m") else "기타", head)
with io.open(OUT + "/tables/BONOBO_READ_BASECALLER.tsv", "w", encoding="utf8", newline="\n") as f:
    f.write("read_id\tbasecaller\tmm_head\n")
    for k in sorted(rows):
        f.write("%s\t%s\t%s\n" % (k, rows[k][0], rows[k][1]))
print(collections.Counter(v[0] for v in rows.values()))
