# -*- coding: utf-8 -*-
"""A1: move the ft nucleosome/MSP calls onto read-unit tracks (unit_pos 1..unit_len).

Input   ft_pilot\\<sample>.<cfg>.<src>.nuc.bam   (ft add-nucleosomes v0.13 defaults)
Convention (measured in p03)
  Ma:Z:<len>;nuc.:<s1>-<l>,...;msp.:... : starts are 1-based, in the original read orientation (as MM)
  on a reverse alignment the SEQ coordinate is sp = len-1-op; SEQ to canonical by CIGAR (M/=/X); canonical to unit by MAP_V11
Units are only the (read, slot, array_pos) with in_final=Y in M6A_READUNIT.
Output  a01_tracks\\<sample>.<src>.tsv.gz
  sample src slot array_pos read_id unit_len track
  track, one character per unit_pos 1.. : N nucleosome, M MSP, . no call (45 bp at the ends, deletion, outside the read)
usage: python a01_unit_tracks.py <path to nuc.bam> ...
"""
import io, os, sys, gzip, bisect, collections
# The data and output roots are environment variables with no default.
#  Paths to files inside this repository are resolved relative to `__file__` - `a_lib`, `bamlite` and `p02_lib` are deposited in this same folder.
#  The table of variables is in `README.md`, section 'Environment variables'.
_HERE = os.path.dirname(os.path.abspath(__file__))


def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s. "
            "The deposited copy has no default; see README.md, section 'Environment variables'."
            % (k, os.path.basename(__file__)))

sys.path.insert(0, _HERE)
from bamlite import read_bam
from p02_lib import seq2ref

R = _need("RDNA_V12_CANON_DIR")
#  This is the folder it writes to - pointing it at a published input folder would overwrite those tables.
OUT = _need("RDNA_A01_TRACKS_OUT_DIR")
os.makedirs(OUT, exist_ok=True)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8", errors="replace",
                              line_buffering=True)
#  File-name labels (RFC1918 private range; not addresses) - can be overridden with `NODE_LABELS`
HOST = {"w3": "ip-172-31-58-58", "iB": "ip-172-31-32-12", "iA": "ip-172-31-57-143"}
if os.environ.get("NODE_LABELS"):
    HOST = dict(kv.split("=", 1) for kv in os.environ["NODE_LABELS"].split(",") if "=" in kv)
METH = {"UCSC_ONT_FIBERSEQ": "fsmeth_v12c", "PACBIO_FIBERSEQ": "pbmeth_v12d"}

UNITS = collections.defaultdict(list)
with io.open(os.path.join(R, "02_canonical", "MAP_V11.tsv"), encoding="utf8") as f:
    ix = {c: i for i, c in enumerate(f.readline().rstrip("\n").split("\t"))}
    for l in f:
        x = l.rstrip("\n").split("\t")
        if x[ix["seg_type"]] == "UNIT":
            UNITS[(x[ix["sample"]], x[ix["slot"]])].append(
                (int(x[ix["unit_start"]]), int(x[ix["unit_end"]]), x[ix["array_pos"]]))


_FINAL = None


def final_units(sample, src):
    """M6A_READUNIT in_final=Y (read, slot, array_pos), read once and kept."""
    global _FINAL
    if _FINAL is None:
        _FINAL = collections.defaultdict(set)
        for s_, meth in METH.items():
            for node in ("w3", "iB", "iA"):
                p = os.path.join(R, "08_methylation_s5", node, meth, "M6A_READUNIT_%s.tsv.gz" % HOST[node])
                with gzip.open(p, "rt", encoding="utf8") as f:
                    f.readline()
                    for l in f:
                        x = l.split("\t", 13)                   # 12 = in_final
                        if x[1] == s_ and x[12] == "Y":
                            _FINAL[(x[0], s_)].add((x[6], x[4], x[5]))
    return _FINAL.get((sample, src), set())


def run(bam):
    base = os.path.basename(bam)
    sample, cfg, src = base.split(".")[:3]
    want = final_units(sample, src)
    c = collections.Counter()
    out = gzip.open(os.path.join(OUT, "%s.%s.tsv.gz" % (sample, src)), "wt", encoding="utf8")
    out.write("sample\tsrc\tslot\tarray_pos\tread_id\tunit_len\ttrack\n")
    done = set()
    for r in read_bam(bam):
        if r["flag"] & 0x900:
            continue
        sec = {}
        for part in r["tags"].get("Ma", "").split(";")[1:]:
            if ":" in part:
                nm, lst = part.split(":", 1)
                sec[nm.rstrip(".+-")] = [(int(a) - 1, int(b)) for a, b in
                                         (v.split("-") for v in lst.split(",") if v)]
        if not sec:
            c["read_no_calls"] += 1
            continue
        n = len(r["seq"]); rev = bool(r["flag"] & 16)
        lab = bytearray(b"." * n)                   # SEQ coordinates
        for key, ch in (("nuc", 78), ("msp", 77)):
            for a, l in sec.get(key, []):
                for op in range(a, min(a + l, n)):
                    lab[(n - 1 - op) if rev else op] = ch
        smp, slot = r["ref"].rsplit("_h", 1)
        us = sorted(UNITS[(smp, slot)])
        ustart = [u[0] for u in us]
        tracks = {}
        for sp, g in seq2ref(r).items():
            j = bisect.bisect_right(ustart, g) - 1
            if j >= 0:
                a, b, ap = us[j]
                if a <= g <= b and (r["name"], slot, ap) in want:
                    t = tracks.setdefault(ap, bytearray(b"." * (b - a + 1)))
                    t[g - a] = lab[sp]
        for ap, t in tracks.items():
            k = (r["name"], slot, ap)
            if k in done:
                c["dup_unit"] += 1
                continue
            done.add(k)
            out.write("%s\t%s\t%s\t%s\t%s\t%d\t%s\n" % (sample, src, slot, ap, r["name"], len(t),
                                                         t.decode()))
            c["unit_written"] += 1
            c["bp"] += len(t); c["N"] += t.count(b"N"); c["M"] += t.count(b"M")
    out.close()
    print("%s %s  final units %d, tracks %d, missing %d, reads without calls %d, duplicates %d, N %.3f M %.3f none %.3f"
          % (sample, src, len(want), c["unit_written"], len(want - done), c["read_no_calls"],
             c["dup_unit"], c["N"] / max(1., c["bp"]), c["M"] / max(1., c["bp"]),
             1 - (c["N"] + c["M"]) / max(1., c["bp"])))


if __name__ == "__main__":
    for b in sys.argv[1:]:
        run(b)
