# -*- coding: utf-8 -*-
"""Shared loader for the a0x scripts, with the definitions of a02_tss_profile.py (canonical TSS, β, track to canonical array)."""
import io, os, gzip, glob, collections
import numpy as np

# The data and output roots are environment variables with no default (CONVENTIONS, "the deposited copy wins").
#  The table of variables is in `README.md`, section 'Environment variables'.

def _need(k):
    try:
        return os.environ[k]
    except KeyError:
        raise SystemExit(
            "environment variable %s is not set; it is required by %s "
            "(a module stops the script that imports it). "
            "The deposited copy has no default; see README.md, section 'Environment variables'."
            % (k, os.path.basename(__file__)))

R = _need("RDNA_V12_CANON_DIR")
D = _need("RDNA_ACCESS_DIR")
#  These three are not addresses but machine labels that appear in the names of the deposited data files
#  (`METH_READUNIT_<label>.tsv.gz`). `172.31.x.x` is the RFC1918 private range and is not reachable.
#  To use data deposited under other names, override them with `NODE_LABELS` (`w3=...,iB=...,iA=...`).
HOST = {"w3": "ip-172-31-58-58", "iB": "ip-172-31-32-12", "iA": "ip-172-31-57-143"}
if os.environ.get("NODE_LABELS"):
    HOST = dict(kv.split("=", 1) for kv in os.environ["NODE_LABELS"].split(",") if "=" in kv)
METH = {"UCSC_ONT_FIBERSEQ": "fsmeth_v12c", "PACBIO_FIBERSEQ": "pbmeth_v12d"}
COMP = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
PB_EXCL = {"HG02630", "hg002"}          # PacBio fibertools 0.4.2 - excluded from the jasmine comparison
BCUT = 0.5


STATE_MODE = os.environ.get("STATE_MODE", "")            # "" = beta < 0.5 - "gmm" = the a13 mixture-model call
SUF = "_gmm" if STATE_MODE == "gmm" else ""


def load_state_gmm():
    """(sample, src, read, slot, apos) -> 'L'/'H'/'U'  (a13_unit_state_gmm, P(L)>=0.9 / <=0.1)"""
    S = {}
    with gzip.open(os.path.join(D, "a13_unit_state_gmm.tsv.gz"), "rt", encoding="utf8") as f:
        f.readline()
        for l in f:
            x = l.rstrip("\n").split("\t")
            S[(x[1], x[0], x[2], x[3], x[4])] = x[10]
    return S


def load_units():
    """(sample, slot, apos) -> (unit_start, unit_end, tss; canonical, 1-based), and the dict of canonical sequences"""
    SUM = {}
    for l in io.open(os.path.join(R, "02_canonical", "SUM_V11.tsv"), encoding="utf8").readlines()[1:]:
        x = l.rstrip("\n").split("\t"); SUM[(x[0], x[1])] = x[4]
    FN = {}
    for l in io.open(os.path.join(R, "02_canonical", "hap_seq_node_v9.tsv"), encoding="utf8"):
        x = l.rstrip("\n").split("\t")
        if len(x) >= 3:
            FN[(x[0], x[1])] = x[2]
    UNIT, SEQ = {}, {}
    with io.open(os.path.join(R, "02_canonical", "MAP_V11.tsv"), encoding="utf8") as f:
        ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
        for l in f:
            x = l.rstrip("\n").split("\t")
            if x[ix["seg_type"]] != "UNIT":
                continue
            k = (x[ix["sample"]], x[ix["slot"]])
            if k not in SEQ:
                s = "".join(z.strip() for z in io.open(os.path.join(R, "02_canonical", "aog8v12fs", "cmp4d_v7", FN[k]))
                            if not z.startswith(">")).upper()
                SEQ[k] = "".join(COMP.get(b, "N") for b in reversed(s)) if SUM.get(k) == "Y" else s
            us, ue = int(x[ix["unit_start"]]), int(x[ix["unit_end"]])
            u = SEQ[k][us - 1:ue]
            p = u.rfind("GTCTACGGC")
            if p < 0:
                p = len(u) - 119
            UNIT[(k[0], k[1], x[ix["array_pos"]])] = (us, ue, us + p)
    return UNIT, SEQ


def load_beta(min_cpg=10):
    """(sample, src, read, slot, apos) -> β  (in_final=Y)"""
    B = {}
    for src, meth in METH.items():
        for node in ("w3", "iB", "iA"):
            p = os.path.join(R, "08_methylation_s5", node, meth, "METH_READUNIT_%s.tsv.gz" % HOST[node])
            with gzip.open(p, "rt", encoding="utf8") as f:
                ix = {k: i for i, k in enumerate(f.readline().rstrip("\n").split("\t"))}
                for l in f:
                    x = l.split("\t", 12)
                    if x[ix["data_source"]] != src or x[ix["in_final"]] != "Y":
                        continue
                    if int(x[ix["n_cpg"]]) >= min_cpg:
                        B[(x[ix["sample"]], src, x[ix["read_id"]], x[ix["slot"]], x[ix["array_pos"]])] = float(x[ix["beta"]])
    return B


def load_m6a():
    """(sample, src, read, slot, apos) -> [n_a, n_m6a]  (in_final=Y, both strands summed)"""
    M = collections.defaultdict(lambda: [0, 0])
    for src, meth in METH.items():
        for node in ("w3", "iB", "iA"):
            p = os.path.join(R, "08_methylation_s5", node, meth, "M6A_READUNIT_%s.tsv.gz" % HOST[node])
            with gzip.open(p, "rt", encoding="utf8") as f:
                f.readline()
                for l in f:
                    x = l.split("\t", 13)
                    if x[1] != src or x[12] != "Y":
                        continue
                    v = M[(x[0], src, x[6], x[4], x[5])]
                    v[0] += int(x[9]); v[1] += int(x[10])
    return M


def iter_reads(UNIT, pad_left, pad_right):
    """per read: (smp, src, rid, slot, arr (uint8 canonical array), lo, [(apos, T_index)])
    arr[i] = 'M'/'N'/'.' at canonical 0-based position lo+i"""
    for tf in sorted(glob.glob(os.path.join(D, "a01_tracks", "*.tsv.gz"))):
        reads = collections.defaultdict(list)
        with gzip.open(tf, "rt", encoding="utf8") as f:
            f.readline()
            for l in f:
                smp, src, slot, ap, rid, ul, tr = l.rstrip("\n").split("\t")
                reads[(smp, src, rid, slot)].append((ap, tr))
        for (smp, src, rid, slot), us in reads.items():
            lo = min(UNIT[(smp, slot, ap)][0] for ap, _ in us) - 1 - pad_left - 1
            hi = max(UNIT[(smp, slot, ap)][1] for ap, _ in us) + pad_right + 1
            arr = np.full(hi - lo, ord("."), dtype=np.uint8)
            for ap, tr in us:
                s0 = UNIT[(smp, slot, ap)][0] - 1 - lo
                arr[s0:s0 + len(tr)] = np.frombuffer(tr.encode(), dtype=np.uint8)
            yield smp, src, rid, slot, arr, lo, [(ap, UNIT[(smp, slot, ap)][2] - 1 - lo) for ap, _ in us]
