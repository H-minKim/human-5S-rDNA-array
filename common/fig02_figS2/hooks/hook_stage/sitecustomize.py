# -*- coding: utf-8 -*-
r"""Hook for recovering the canonical export. It does not touch the variant reference (CHM13 as is).
   It does one thing: redirect G.OUT to staging so nothing can be written to the canonical output.
   STAGEF (input data) stays canonical, it is only read.
"""
import os, sys
STAGE = os.environ["STAGE_OUT"]
sys.path.insert(0, __import__("os").environ.get("V11_CODE", __import__("os").path.join(__import__("os").path.dirname(__import__("os").path.abspath(__file__)), __import__("os").pardir, __import__("os").pardir)))
import v11_io as G
G.OUT = STAGE
assert os.path.normcase(os.path.abspath(G.OUT)) != os.path.normcase(
    os.path.abspath((__import__("os").environ["RDNA_DATA_ROOT"] + r"/5S_position_variant_work/2026-09-17_v11_fix2"))), "정본 경로 방지"
