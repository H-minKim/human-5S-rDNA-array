# -*- coding: utf-8 -*-
# Korean edition: hook_tier (first) + symbol cleanup hook
import os, sys
_H = os.path.dirname(os.path.abspath(__file__))
_T = os.path.join(os.path.dirname(_H), "hook_tier")
sys.path.insert(0, _T)
import importlib.util as _u
_s = _u.spec_from_file_location("_hook_tier_sc", os.path.join(_T, "sitecustomize.py"))
_m = _u.module_from_spec(_s); _s.loader.exec_module(_m)
sys.path.insert(0, _H)
import _symclean
_symclean.install()
