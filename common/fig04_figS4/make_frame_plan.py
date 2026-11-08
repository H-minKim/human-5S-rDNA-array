# -*- coding: utf-8 -*-
# 틀 계획을 낸다 - 그림마다 frame_layout.plan 을 돌려 SPEC_FRAME json 으로.
import io, json, os, re, sys
SP = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, SP)
NL = chr(10)
import frame_layout
from frame_layout import plan
# S4 의 D 는 지금의 80%, E,F 는 더 길게, B 도 더 길게
frame_layout.WEIGHT = {
    'ClaimA_two_states__c': 0.80,
    'Fig4c2c3c4_block_anatomy__c3': 1.30,
    'Fig4c2c3c4_block_anatomy__c4': 1.30,
    'Fig4d_single_molecules__m1': 1.18, 'Fig4d_single_molecules__m2': 1.18,
    'Fig4d_single_molecules__m3': 1.18, 'Fig4d_single_molecules__m4': 1.18,
}

# 배치표는 조립기(v9)에서 그대로 읽는다 - 거기 한 벌만 둔다
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

# A,B,C 를 세로로 늘린다
for _k, _h in (("Fig4_a_a2_b_c2_c3_v3__a", 78.0),
               ("Fig4_a_a2_b_c2_c3_v3__b", 78.0),
               ("Fig4_a_a2_b_c2_c3_v3__a2", 62.0)):
    if _k in main_plan:
        main_plan[_k][1] = _h
# , E,F 는 지금의 70% , D 는 세로로 늘림 , G~O 는 C 와 같은 세로(62)
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

# F(축이 둘이라 준 값의 0.818 배로 그려진다) 는 되돌려 적는다
if "Claim_G_outside_array__b" in main_plan:
    main_plan["Claim_G_outside_array__b"][1] = round(
        main_plan["Claim_G_outside_array__b"][1] / 0.818, 1)

# S6 DXZ4 줄 - H 를 I 와 같은 가로로,
# G 가 남은 가로를 차지한다
_h, _i, _g = ("F39_DXZ4_array__b", "F39_DXZ4_array__b-2", "F39_DXZ4_array__c")
if all(k in main_plan for k in (_h, _i, _g)):
    _d = main_plan[_h][0] - main_plan[_i][0]
    main_plan[_h][0] = main_plan[_i][0]
    main_plan[_g][0] = round(main_plan[_g][0] + _d, 3)

io.open(os.path.join(SP, "frame_main.json"), "w", encoding="utf8", newline=NL).write(
    json.dumps(main_plan, ensure_ascii=False, indent=1))
# S4 의 세로축 이름이 틀보다 길어 패널 글자와 부딪혔다
# (Fitted component mean β 55 pt ↔ 틀 44.2) ⇒ 글을 깎지 말고 틀을 키운다
for _k, _h in (("SuppS6_from_c", 56.0), ("ClaimA_two_states__b", 56.0),
               ("Fig4d_single_molecules__m1", 56.0), ("Fig4d_single_molecules__m2", 56.0),
               ("Fig4d_single_molecules__m3", 56.0), ("Fig4d_single_molecules__m4", 56.0),
               ("ClaimA_two_states__c", 58.0),
               ("Fig4c2c3c4_block_anatomy__c3", 58.0),
               ("Fig4c2c3c4_block_anatomy__c4", 58.0)):
    if _k in s41_plan:
        s41_plan[_k][1] = _h

# S4 의 K 는 범례가 틀을 벗어나지 않게 가로를 키운다
# (범례 2 칸 실측 폭 102 pt ↔ 옛 틀 80.79). 줄의 틀 합은 그대로 두고 나머지 둘에서 뺀다.
_I, _J, _K = ("F39_5S_array__e", "F39_D4Z4_4q_array__e",
              "Claim_F2_two_modes_of_transition__c")
if all(k in s41_plan for k in (_I, _J, _K)):
    _tot = s41_plan[_I][0] + s41_plan[_J][0] + s41_plan[_K][0]
    s41_plan[_K][0] = 104.0
    s41_plan[_I][0] = s41_plan[_J][0] = round((_tot - 104.0) / 2.0, 3)

io.open(os.path.join(SP, "frame_s41.json"), "w", encoding="utf8", newline=NL).write(
    json.dumps(s41_plan, ensure_ascii=False, indent=1))
print()
print("낸 것 : frame_main.json %d 칸 · frame_s41.json %d 칸"
      % (len(main_plan), len(s41_plan)))
