# -*- coding: utf-8 -*-
# English version: hook_tier -> install symbol cleanup -> install translation
#   The hook installed last runs first. The translation must be installed last so that
#     the order becomes 'match the dictionary on the Korean original -> clean symbols in the English result'.
#     (the dictionary keys contain stars and middle dots, so reversing the order makes all 945 entries fail to match)
import os, sys
_H = os.path.dirname(os.path.abspath(__file__))
_B = os.path.dirname(_H)
_T = os.path.join(_B, "hook_tier")
sys.path.insert(0, _T)
import importlib.util as _u
_s = _u.spec_from_file_location("_hook_tier_sc", os.path.join(_T, "sitecustomize.py"))
_m = _u.module_from_spec(_s); _s.loader.exec_module(_m)

sys.path.insert(0, os.path.join(_B, "hook_sym"))
import _symclean
_symclean.install()

sys.path.insert(0, os.path.join(_B, "hook_en"))
import _entrans
_entrans.install()
if os.environ.get("EN_WIDTH_LOG"):
    import _width
    _width.install()
