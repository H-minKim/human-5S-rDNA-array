# -*- coding: utf-8 -*-
# PERARRAY_META_v3 의 유전자 구간을 strand 보정해서 다시 계산해 고친다
#   (x42 가 유전자 구간을 g8 틀로 뒤집지 않아, strand '-' array 묶음에서 음영이 반대쪽에 찍혔다. BAM 재스캔은 필요 없다.)
#   유닛 안 좌표 : + array → x - unit_start + 1 · − array → unit_end - y + 1  (유전자가 유닛 끝에 오는 g8 틀)
import io, os, sys, glob, re
import numpy as np
import pandas as pd

import os as _os_need, sys as _sys_need


def _need(k):
    """No defaults are set: if an environment variable is missing the script stops instead of silently reading the wrong location."""
    v = _os_need.environ.get(k)
    if not v:
        _sys_need.exit("environment variable %s is not set (the deposited copy has no default; "
                       "see README.md, section 'How to run it')" % k)
    return v


NEW = _need("RDNA_OUT_DIR")
B = _need("RDNA_DATA_ROOT")
OLD = _need("RDNA_PREV_OUT_DIR")
sys.path.insert(0, NEW + "/code")
import arch  # noqa
try:
    sys.stdout.reconfigure(encoding="utf-8")
except Exception:
    pass
src = io.open(OLD + "/code/x11_coding_readunit.py", encoding="utf8").read()
L11 = 'LOG = io.open(OUT + ("/X11_CODING_PILOT.txt" if PILOT else "/X11_CODING_BUILD.txt"), "a", encoding="utf8", newline="\\n")'
head = src.split("STAT = collections.Counter()")[0].replace(L11, "LOG = io.StringIO()")
G = {"__name__": "h"}
exec(compile(head, "x11", "exec"), G)
load_species = G["load_species"]
M = pd.read_csv(NEW + "/tables/PERARRAY_META_v3.tsv", sep="\t")
SPK = {"고릴라": "gorilla_mGorGor1", "보노보": "bonobo_mPanPan1", "침팬지": "chimp_mPanTro3",
       "보르네오오랑우탄": "orang_bornean_mPonPyg2", "수마트라오랑우탄": "orang_sumatran_mPonAbe1",
       "시아망": "siamang_mSymSyn1", "애기장대": "arabidopsis_ColCC",
       "금화조": "zfinch_v03_dip", "제브라피시": "zebrafish_GRCz12tu"}
CACHE, STR = {}, {}
for ko, sp in SPK.items():
    CACHE[sp] = load_species(sp)
    for f in sorted(glob.glob(B + "/02_consensus_map/%s/*.g8map.tsv" % sp)):
        if not re.search(r"_[0-9]+_[0-9]+[.]g8map[.]tsv$", f.replace("\\", "/")):
            continue
        g = pd.read_csv(f, sep="\t", usecols=["array_id", "strand"])
        STR[g.array_id.iloc[0]] = g.strand.iloc[0]
old_s, old_e, new_s, new_e, flipped = [], [], [], [], 0
for t in M.itertuples(index=False):
    sp = SPK[t.종]
    aids = str(t.array목록).split("|")
    gw = []
    for aid in aids:
        v = CACHE[sp][aid]
        st = STR[aid]
        for i in range(len(v["us"])):
            for (x, y) in v["genes"][i]:
                gw.append((x - v["us"][i] + 1, y - v["us"][i] + 1) if st == "+" else (v["ue"][i] - y + 1, v["ue"][i] - x + 1))
    s_ = float(np.median([a for a, b in gw])) if gw else np.nan
    e_ = float(np.median([b for a, b in gw])) if gw else np.nan
    old_s.append(t.유전자시작); old_e.append(t.유전자끝); new_s.append(s_); new_e.append(e_)
    if np.isfinite(s_) and np.isfinite(t.유전자시작) and abs(s_ - t.유전자시작) > 1:
        flipped += 1
M["유전자시작"], M["유전자끝"] = new_s, new_e
TMP = os.environ.get("TEMP", ".") + "/x47"
os.makedirs(TMP, exist_ok=True)
M.to_csv(TMP + "/PERARRAY_META_v3.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/PERARRAY_META_v3.tsv", "tables/PERARRAY_META_v3.tsv", "유전자 구간 strand 보정(2026-10-01)")
D = pd.DataFrame(dict(종=M.종, 묶음=M.묶음, 유닛길이=M.유닛길이, 옛시작=old_s, 옛끝=old_e, 새시작=new_s, 새끝=new_e))
D["바뀜"] = np.where(np.abs(D.새시작 - D.옛시작) > 1, "바뀜", "같음")
D.to_csv(TMP + "/GENE_WINDOW_FIX_v3.tsv", sep="\t", index=False, encoding="utf8", lineterminator="\n")
arch.put(TMP + "/GENE_WINDOW_FIX_v3.tsv", "tables/GENE_WINDOW_FIX_v3.tsv", "유전자 구간 strand 보정 전후(2026-10-01)")
print("묶음 %d 중 %d 개의 유전자 구간이 바뀜 (strand '-' 묶음)" % (len(M), flipped))
print(D[D.바뀜 == "바뀜"].head(12).to_string(index=False))
