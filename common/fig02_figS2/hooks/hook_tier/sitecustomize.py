# -*- coding: utf-8 -*-
# ============================================================================
# hook_tier = hook_majref + cohort tier swap
#   the same as hook_majref except for the one thing below
#   if V11_TIER = h1 | h3, G.load("primary") is replaced by loading that preset.
#   => the tier is applied without editing the 73 scripts that hardcode "primary".
#   If V11_TIER is not set, nothing changes (= identical to hook_majref).
#   Output paths are confined to the mirror by REPRO_MIRROR, so the canonical copy (F:) cannot be written.
# ============================================================================
r"""Replaces v11_io.load so that the majority-form consensus is used as the reference

Why: the current reference `unit_canon.fa` = the revcomp single unit of CHM13 chr1:227,887,698-227,889,937.
     When the majority form was chosen at all 2,240 positions, CHM13 carried a minority form at 11 positions, each with
     one token carried by more than half of the units (11 tokens in all), inflating every similarity measure.
     Using the majority form at those 11 positions removes 8.64 (55.7%) of the 15.515 variants per unit in the
     564-haplotype set, and 8.65 (55.8%) of the 15.520 in H1 (557 haplotypes).
     Position 1476 has a most frequent form (IB46:2, 0.42 of units, 29 alleles) but no majority, so it is not flipped.

What changes: only the variant set of a unit (the 5th item of a row).
   position p not flipped      : tokens unchanged
   position p flipped (11)     : a variant only when the unit allele differs from M(p)
        - a different token t -> t          - no token (= CHM13 form) -> "REF@p"
meta['raw'] is not touched (used only by EXPORT_variant_groups.py, do not use that script with this version)
Canonical path protection: if REPRO_MIRROR is set, G.OUT is also redirected to the mirror
"""
import os, sys, re, collections

MIRROR = os.environ.get("REPRO_MIRROR")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir, __import__("os").pardir)))
try:
    import v11_io as G
except Exception as e:
    sys.stderr.write("sitecustomize: v11_io import failed %r\n" % e); raise

if MIRROR:
    G.OUT = MIRROR
    G.STAGEF = os.path.join(MIRROR, "STAGE_ALL_V11F2.tsv")
    assert not G.OUT.startswith(("F:", "f:")), "canonical path guard"

BLOCKS = G.BLOCKS
def _pos(t):
    if t.startswith("DB") or t.startswith("IB"):
        b = int(t[2:].split(":")[0]); return BLOCKS[b - 1][0] if 1 <= b <= len(BLOCKS) else None
    if t[0] == "D": return int(t[1:].split(":")[0])
    if t[0] == "I": return int(t[1:].split(":")[0])
    if t[0] == "R": return int(t[1:-2])
    try: return int(t[:-2])
    except Exception: return None

_ORIG = G.load
_CACHE = {}

def _majority(U):
    alle = collections.defaultdict(collections.Counter); nu = 0
    for k in U:
        for r in U[k]:
            nu += 1; dd = {}
            for t in r[4]:
                p = _pos(t)
                if p is None: continue
                dd[p] = min(dd[p], t) if p in dd else t
            for p, t in dd.items(): alle[p][t] += 1
    for p in alle: alle[p]["REF"] = nu - sum(alle[p].values())
    FLIP = {}
    for p in alle:
        t, c = alle[p].most_common(1)[0]
        # change only when the majority form is not REF and exceeds half of the units
        if t != "REF" and c / nu > 0.5: FLIP[p] = t
    return FLIP, nu

def _inv(t):
    """Inverse of the majority tokens, the variants the 'CHM13 form' has against the majority reference"""
    m = re.match(r"^(R?)(\d+)([ACGT])([ACGT])$", t)
    if m: return "%s%s%s%s" % (m.group(1), m.group(2), m.group(4), m.group(3))
    m = re.match(r"^D(\d+):(\d+)$", t)
    if m: return "I%s:%s:?" % (m.group(1), m.group(2))
    m = re.match(r"^I(\d+):(\d+):.*$", t)
    if m: return "D%s:%s" % (m.group(1), m.group(2))
    m = re.match(r"^DB(\d+):(\d+)$", t)
    if m: return "IB%s:%s" % (m.group(1), m.group(2))
    m = re.match(r"^IB(\d+):(\d+)$", t)
    if m: return "DB%s:%s" % (m.group(1), m.group(2))
    raise ValueError("cannot build the reverse token: %r" % t)

_TIER = os.environ.get("V11_TIER", "").strip()
if _TIER:
    sys.stderr.write("[TIER] V11_TIER=%s — G.load('primary') is loaded as this preset\n" % _TIER)
    assert _TIER in G.PRESETS, "unknown preset %r" % _TIER
    assert MIRROR, "REPRO_MIRROR must be set when V11_TIER is used (canonical path protection)"

def load(preset="primary", *a, **kw):
    if _TIER and preset == "primary":
        preset = _TIER
    U, A, M = _ORIG(preset, *a, **kw)
    key = preset
    if key not in _CACHE:
        _CACHE[key] = _majority(U)
    FLIP, nu = _CACHE[key]
    n_old = n_new = 0
    for k in list(U):
        rows = []
        for r in U[k]:
            dd = {}
            for t in r[4]:
                p = _pos(t)
                if p is None: continue
                dd[p] = min(dd[p], t) if p in dd else t
            s = set()
            for p, t in dd.items():
                if p not in FLIP or t != FLIP[p]: s.add(t)
            for p, m in FLIP.items():
                if p not in dd: s.add(_inv(m))
            n_old += len(r[4]); n_new += len(s)
            rows.append(tuple(r[:4]) + (s,) + tuple(r[5:]))
        U[k] = rows
    sys.stderr.write("[MAJREF] preset=%s · flipped positions %d · variants per unit %.3f -> %.3f\n"
                     % (preset, len(FLIP), n_old / max(nu, 1), n_new / max(nu, 1)))
    return U, A, M

G.load = load

# ---- Fig 2 and S2 on one set: keep only the (sample, slot) pairs that V11_SUBSET points to ----
# The majority reference (FLIP) is set from the cohort before filtering, so the reference is unchanged. Without the environment variable nothing happens.
_SUBF = os.environ.get("V11_SUBSET")
if _SUBF:
    _KEEP = set()
    for _l in open(_SUBF, encoding="utf8"):
        _a = _l.rstrip("\n").split("\t")
        if len(_a) == 2: _KEEP.add((_a[0], _a[1]))
    _LOAD0 = G.load
    def _load_sub(*a, **kw):
        U, A, M = _LOAD0(*a, **kw)
        U2 = {k: v for k, v in U.items() if (str(k[0]), str(k[1])) in _KEEP}
        if len(U2) != len(_KEEP):
            raise SystemExit("V11_SUBSET is not inside the load result %d %d" % (len(U2), len(_KEEP)))
        sys.stderr.write("[SUBSET] %d -> %d\n" % (len(U), len(U2)))
        return U2, A, M
    G.load = _load_sub
