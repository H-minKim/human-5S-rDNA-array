# -*- coding: utf-8 -*-
# ============================================================================
# hook_tieronly = cohort tier swap only. The majority (majref) conversion is not done.
#
#   Why: MAJREF.py is code that computes the majority form itself.
#        with hook_tier, G.load already returns data converted to the majority form, so
#        MAJREF.py would apply it twice (measured: 57% reduction at 12 positions -> 3% reduction at 1 position).
#   So fig_T0_majref and fig_unitpos_majref are run with this hook.
#   Canonical protection: if REPRO_MIRROR is set, G.OUT is redirected to the mirror (same as hook_tier).
# ============================================================================
import os, sys

MIRROR = os.environ.get("REPRO_MIRROR")
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir, __import__("os").pardir)))
import v11_io as G

if MIRROR:
    G.OUT = MIRROR
    G.STAGEF = os.path.join(MIRROR, "STAGE_ALL_V11F2.tsv")
    assert not G.OUT.startswith(("F:", "f:")), "정본 경로 방지"

_TIER = os.environ.get("V11_TIER", "").strip()
_ORIG = G.load
if _TIER:
    sys.stderr.write("[TIERONLY] V11_TIER=%s — preset 만 바꾼다 (다수형 변환 안 함)\n" % _TIER)
    assert _TIER in G.PRESETS, "모르는 preset %r" % _TIER
    assert MIRROR, "V11_TIER 를 쓸 때는 REPRO_MIRROR 가 있어야 한다"

    def load(preset="primary", *a, **kw):
        if preset == "primary":
            preset = _TIER
        return _ORIG(preset, *a, **kw)

    G.load = load

# symbol cleanup hook
_H = os.path.dirname(os.path.abspath(__file__))          # .../scratchpad/hook_tieronly
sys.path.insert(0, os.path.join(os.path.dirname(_H), "hook_sym"))
import _symclean
_symclean.install()
