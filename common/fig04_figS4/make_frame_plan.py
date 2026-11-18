# -*- coding: utf-8 -*-
# Produce the frame plan: frame_layout.plan is run per figure and written as a SPEC_FRAME json.
import io, json, os, re, sys
SP = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SP)
NL = chr(10)
import frame_layout
from frame_layout import plan
# horizontal weights: D of S4 is narrowed to 0.80, and E, F and B are lengthened
frame_layout.WEIGHT = {
    'ClaimA_two_states__c': 0.80,
    'Fig4c2c3c4_block_anatomy__c3': 1.30,
    'Fig4c2c3c4_block_anatomy__c4': 1.30,
    'Fig4d_single_molecules__m1': 1.18, 'Fig4d_single_molecules__m2': 1.18,
    'Fig4d_single_molecules__m3': 1.18, 'Fig4d_single_molecules__m4': 1.18,
}

# the layout table is read from the assembler (version 9), where the only copy is kept
src = io.open(os.path.join(SP, "assemble_v9_2026-10-18.py"), encoding="utf8").read()
ns = {"__file__": os.path.join(SP, "assemble_v9_2026-10-18.py")}
head = src[:src.index("def fetch(")]
exec(compile(head, "assemble_head", "exec"), ns)
FIG4, S41, S42, S44 = ns["FIG4"], ns["S41"], ns["S42"], ns["S44"]

SPEC = os.path.join(SP, "panelsplit_spec", "_panels")
S41D = os.path.join(SP, "panelsplit_s41", "_panels")

main_plan = {}
main_plan.update(plan(FIG4, SPEC, "Fig4_main"))
main_plan.update(plan(S44, SPEC, "SuppFig_S6"))
s41_plan = plan(S41, S41D, "SuppFig_S4")

# A, B and C are stretched vertically
for _k, _h in (("Fig4_a_a2_b_c2_c3_v3__a", 78.0),
               ("Fig4_a_a2_b_c2_c3_v3__b", 78.0),
               ("Fig4_a_a2_b_c2_c3_v3__a2", 62.0)):
    if _k in main_plan:
        main_plan[_k][1] = _h
# E and F are reduced in height, D is stretched vertically, and G to O take the same height as C, 62
for _k, _h in (("Fig4c2c3c4_block_anatomy__c2", 39.4),
               ("Claim_G_outside_array__b", 39.4),
               ("Fig4c_six_units_continuous__c", 68.0),
               ("F39_5S_array__a", 62.0), ("F39_5S_array__c", 62.0),
               ("F39_5S_array__b", 62.0),
               ("F39_D4Z4_4q_array__a", 62.0), ("F39_D4Z4_4q_array__c", 62.0),
               ("F39_D4Z4_4q_array__b", 62.0),
               ("Claim_F2_show_hardcall_ext_EN__d10", 62.0),
               ("Claim_F2_show_hardcall_ext_EN__d9", 62.0),
               ("Claim_F2_show_hardcall_ext_EN__d11", 62.0)):
    if _k in main_plan:
        main_plan[_k][1] = _h

# F is drawn at 0.818 of the value given, having two axes, so the value is adjusted back
if "Claim_G_outside_array__b" in main_plan:
    main_plan["Claim_G_outside_array__b"][1] = round(
        main_plan["Claim_G_outside_array__b"][1] / 0.818, 1)

# the DXZ4 row of S6: H takes the same width as I,
# and G takes the rest of the width
_h, _i, _g = ("F39_DXZ4_array__b", "F39_DXZ4_array__b-2", "F39_DXZ4_array__c")
if all(k in main_plan for k in (_h, _i, _g)):
    _d = main_plan[_h][0] - main_plan[_i][0]
    main_plan[_h][0] = main_plan[_i][0]
    main_plan[_g][0] = round(main_plan[_g][0] + _d, 3)

io.open(os.path.join(SP, "frame_main.json"), "w", encoding="utf8", newline=NL).write(
    json.dumps(main_plan, ensure_ascii=False, indent=1))
# the y-axis label of S4 is longer than the planned frame and would collide with the panel letter
# ("Fitted component mean beta" at 55 pt against a frame of 44.2), so the frame is widened rather than the text cut
for _k, _h in (("SuppS6_from_c", 56.0), ("ClaimA_two_states__b", 56.0),
               ("Fig4d_single_molecules__m1", 56.0), ("Fig4d_single_molecules__m2", 56.0),
               ("Fig4d_single_molecules__m3", 56.0), ("Fig4d_single_molecules__m4", 56.0),
               ("ClaimA_two_states__c", 58.0),
               ("Fig4c2c3c4_block_anatomy__c3", 58.0),
               ("Fig4c2c3c4_block_anatomy__c4", 58.0)):
    if _k in s41_plan:
        s41_plan[_k][1] = _h

# K of S4 has its width increased so that the legend does not leave the frame
# (the two-column legend measures 102 pt against a planned frame of 80.79). The row's total frame width is unchanged and the difference taken from the other two.
_I, _J, _K = ("F39_5S_array__e", "F39_D4Z4_4q_array__e",
              "Claim_F2_two_modes_of_transition__c")
if all(k in s41_plan for k in (_I, _J, _K)):
    _tot = s41_plan[_I][0] + s41_plan[_J][0] + s41_plan[_K][0]
    s41_plan[_K][0] = 104.0
    s41_plan[_I][0] = s41_plan[_J][0] = round((_tot - 104.0) / 2.0, 3)

io.open(os.path.join(SP, "frame_s41.json"), "w", encoding="utf8", newline=NL).write(
    json.dumps(s41_plan, ensure_ascii=False, indent=1))
print()
print("written: frame_main.json %d cells, frame_s41.json %d cells"
      % (len(main_plan), len(s41_plan)))
