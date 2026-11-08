# -*- coding: utf-8 -*-
# English edition + tier replacement only (no plural conversion): hook_tieronly -> install translation
#   hook_tieronly already installs the symbol cleanup, so the translation is installed later
#   to get the order 'translation -> symbol cleanup'.
import os, sys
_H = os.path.dirname(os.path.abspath(__file__))
_B = os.path.dirname(_H)
_T = os.path.join(_B, "hook_tieronly")
import importlib.util as _u
_s = _u.spec_from_file_location("_hook_tieronly_sc", os.path.join(_T, "sitecustomize.py"))
_m = _u.module_from_spec(_s); _s.loader.exec_module(_m)

sys.path.insert(0, os.path.join(_B, "hook_en"))
import _entrans
_entrans.install()
