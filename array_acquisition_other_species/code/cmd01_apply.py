# -*- coding: utf-8 -*-
# 정정 적용 — 옛 파일은 덮어쓰지 않는다(새 판 *_2026-09-21 을 만든다).
#   ① 극성표 CMD01_POLARITY.tsv  ② 인용 기본값 CMD01_CITE_DEFAULTS.tsv
#   ③ ADULT_DIFFS_2026-09-21.tsv(인용 가능 표시)  ④ README 새 판 6 개(정정 블록 + 옛 문장 표시)
# MIT License · Copyright (c) 2026 Hungmin Kim
import csv
import collections
import glob
import io
import os

R = "$RDNA_OTHERSP_DIR/data_v2/"
T = "\t"
NL = "\n"
V = "2026-09-21"


def wr(p, s):
    assert not os.path.exists(p) or p.endswith(".tsv") and "CMD01_" in p, "덮어쓰기 금지 : " + p
    assert not any(ord(c) < 32 and c not in "\n\t" for c in s), "제어문자"
    io.open(p, "w", encoding="utf8", newline=NL).write(s)


# ① 극성표 ─────────────────────────────────────────────────────────
rows = []
for p in sorted(glob.glob(R + "06_units_from_canonical/*.units.tsv")):
    rr = list(csv.DictReader(io.open(p, encoding="utf8"), delimiter=T))
    if not rr:
        continue
    rr.sort(key=lambda r: int(r["array_pos"]))
    a1, aN = rr[0], rr[-1]
    strand = collections.Counter(r["strand"] for r in rr).most_common(1)[0][0]
    s1, sN = int(a1["chrom_start"]), int(aN["chrom_start"])
    side = "낮은 좌표" if s1 < sN else ("높은 좌표" if s1 > sN else "유닛 1")
    rule_ok = (strand == "+" and side == "낮은 좌표") or (strand == "-" and side == "높은 좌표") or side == "유닛 1"
    hum = ""
    if a1["species_key"] == "human_CHM13v2":
        hum = "array_pos 1 = chr1 높은 좌표 끝 = v11 '게놈 원위(DJ)' 와 같은 끝(v11 좌표 chr1:228,025,900 은 FLANK_CONTEXT 판 · 판 두 벌 있음)"
    rows.append([a1["species_key"], a1["array_id"], a1["chrom"], strand, side, a1["chrom_start"] + "-" + a1["chrom_end"],
                 aN["chrom_start"] + "-" + aN["chrom_end"], str(len(rr)), "Y" if rule_ok else "N", hum])
s = T.join(["species", "array_id", "chrom", "gene_strand", "array_pos1_side", "array_pos1_coords", "array_posN_coords",
            "n_units", "rule_array_pos1_is_transcription_upstream", "human_v11_note"]) + NL
s += NL.join(T.join(r) for r in rows) + NL
wr(R + "00_catalog/CMD01_POLARITY.tsv", s)
print("① 극성표 %d array · 규칙 일치 %d" % (len(rows), sum(1 for r in rows if r[8] == "Y")))

# ② 인용 기본값 ─────────────────────────────────────────────────────
D = [
    ("보노보 CG 메틸", "0.803", "guppy 6.3.8 flowcell read 만 · 호출 가중(Σn_5mC/Σn_calls) · G8 read×유닛", "CMD01_BONOBO_BY_BASECALLER.tsv(재실측)",
     "0.754 = 6.3.7+6.3.8 혼합(PLACEMENT_METH_REPORT) · 0.751 = 09-12 초기 요약(혼합) · 0.736 = AVF README array 중앙(다른 집계) · 0.799 = 옛 문서값(판본별 원표 없음 → 재실측 0.803 으로 교체) · 6.3.7 만 0.667", "종간 비교에는 6.3.8 인 종만"),
    ("침팬지 CG 메틸", "0.646", "전 flowcell(6.3.8) 호출 가중 · G8 read×유닛", "PLACEMENT_METH_REPORT.tsv", "0.641 = 09-12 초기 요약(cross_species_summary)", ""),
    ("시아망 07 합본 일치도", "0.094", "엄격판(치환 + 3 bp 초과 indel) concordance_strict", "07_consensus_from_reads/CONSENSUS_SUMMARY.tsv",
     "0.088 = 비엄격판(옛 인계문 인용은 오기). 다른 종 엄격판 범위 0.696–0.982(제브라 0.696–0.717 · 금화조 UL 0.982)", "07 합본은 엄격판이 규칙"),
    ("금화조 mapq0", "정의 둘", "① chr2 array 별(UL) 89.9–97.0% · TIER_STATUS_FULL 중앙 89.9% = array×시료 중앙 ② 배치 read 전체 UL 69.3% · LIG 71.6%", "TIER_STATUS_FULL.tsv · PLACEMENT_METH_README.md",
     "옛 인계문 '90–97%' 는 ① 의 chr2 만", "분모(chr2 array vs 전체 read)를 반드시 같이 적는다"),
    ("제브라 array 수", "31 (배치·메틸)", "chr4 30 + chr18 1 = read 배치·메틸·AVF 를 잰 array", "CMD01_ARRAY_COUNT_DEFS.tsv",
     "46 = cmsearch array 후보 전체(카탈로그·구조·정본) · 32 = 참조유닛 행(31 + 56528778 의 .k2) · 21 = 참조유닛 '한 종류' 판정 통과 · chr4 6 = WGBS 구역에서 쓴 큰 chr4 array 6 개", "chr4 는 30 array(옛 '31' 은 chr18 포함 총수를 chr4 로 적은 오기)"),
    ("금화조 array 수", "8", "4 좌위 × 2 하플로타입(배치·메틸)", "CMD01_ARRAY_COUNT_DEFS.tsv",
     "33 = cmsearch 후보 전체 · 6 = 서로 다른 서열(chr2 두 쌍이 바이트 동일) · 10 = AVF README 표값(출처 불명 — 원표 AVF.tsv 는 8 array×2 라이브러리)", ""),
    ("애기장대 array 수", "9", "배치·메틸 array", "CMD01_ARRAY_COUNT_DEFS.tsv", "10 = 카탈로그·정본(작은 array 포함) · 5 = 유닛 종류 판정 통과", "CP116284.1 = chr5(README_FIX 의 'chr4' 는 오기)"),
    ("붉은털 array 수", "4", "TIER · 지도(g8map 4)", "CMD01_ARRAY_COUNT_DEFS.tsv", "5 = 카탈로그·정본(NON_5S 1 포함) · 참조유닛 표 4 행 = 3 + NON_5S · NO_UNIT `CM111659.2_103988607` 은 참조유닛 표에 없음", ""),
    ("생쥐 reads_placed", "0 (정의상)", "TIER_STATUS_FULL 의 reads_placed 는 G8 read×유닛 표를 지난 read 수 — 생쥐는 영역 추출 BAM 만 있고 그 표가 없다", "TIER_STATUS_FULL.tsv", "영역 추출 read 127 · 깊이 9.2 · junction 30 은 영역 BAM 에서", "'배치 실패' 가 아니라 '유닛 표 미작성'"),
    ("제브라 배치 read", "3,963–9,267 / 개체", "G8 read×유닛 표를 지난 read", "PLACEMENT_METH_REPORT.tsv", "5,119–11,670 = 09-12 초기 요약 · mapq0 3.5% = array×시료 중앙 vs 11–13% = 개체 전체 read", ""),
    ("애기장대 read", "6,852", "G8 read×유닛 표를 지난 read", "PLACEMENT_METH_REPORT.tsv / TIER", "7,247 재콜링 · 6,972 정렬", ""),
    ("제브라 chr18 유닛·유전자", "유닛 62(OK 59 · NO_FULL 3) · 유전자 65", "정본 유닛표 · cmsearch 유전자 히트(array 안)", "06_units_from_canonical · zfgt.gene_hits.tsv", "64 = 09-12 초기 요약(출처 불명)", ""),
    ("GSC 암컷 코딩−바깥 구간", "−0.076 [−0.163, +0.012]", "GSC_STATS.tsv(5,000 회)", "gsc_GSE122695/GSC_STATS.tsv", "ADULT_DIFFS 는 2,000 회 [−0.164, +0.011] — 같은 결론", "인용은 GSC_STATS"),
]
s = T.join(["item", "default", "definition", "source", "other_values_and_why", "note"]) + NL
s += NL.join(T.join(x) for x in D) + NL
wr(R + "00_catalog/CMD01_CITE_DEFAULTS.tsv", s)
print("② 기본값 %d" % len(D))

# ③ ADULT_DIFFS 새 판 ────────────────────────────────────────────────
st = {}
for r in csv.DictReader(io.open(R + "12_zebrafish_devseries/ADULT_STATS.tsv", encoding="utf8"), delimiter=T):
    st[(r["sample"], r["branch"], r["part"])] = int(r["n_read"])
out = []
hdr = None
for r in csv.DictReader(io.open(R + "12_zebrafish_devseries/ADULT_DIFFS.tsv", encoding="utf8"), delimiter=T):
    hdr = list(r.keys())
    part = r["contrast"].split("-")[0]
    n = st.get((r["sample"], r["branch"], part), 0)
    w = abs(float(r["hi95"]) - float(r["lo95"]))
    if n < 10:
        flag = "인용금지(read<10)"
    elif n < 30:
        flag = "약함(read 10–29)"
    else:
        flag = "가능"
    if n < 10 and w < 0.05:
        flag += " · 구간 좁음은 read 적어 생긴 인공"
    out.append([r[k] for k in hdr] + [str(n), flag])
s = T.join(hdr + ["n_read_part", "cite"]) + NL + NL.join(T.join(x) for x in out) + NL
p = R + "12_zebrafish_devseries/ADULT_DIFFS_" + V + ".tsv"
wr(p, s)
print("③ ADULT_DIFFS 새 판 · 인용금지 %d · 약함 %d" % (sum(1 for x in out if x[-1].startswith("인용금지")), sum(1 for x in out if x[-1].startswith("약함"))))

# ④ README 새 판 ─────────────────────────────────────────────────────
def newver(src, blocks, marks, repl=()):
    s = io.open(src, encoding="utf8").read()
    lines = s.split(NL)
    for key, tag in marks:
        hit = 0
        for i, L in enumerate(lines):
            if key in L and "⟦옛 판" not in L:
                lines[i] = L + "  ⟦옛 판 — " + tag + "⟧"
                hit += 1
        assert hit, "표시할 문장을 못 찾음 : %s · %s" % (src, key)
    s = NL.join(lines)
    for a, b in repl:
        assert a in s, "바꿀 문장 없음 : " + a
        s = s.replace(a, b)
    head = ("> **새 판 " + V + "(정정)** — 옛 판 `" + os.path.basename(src) + "` 는 그대로 둔다(이력용).\n"
            "> 아래 정정이 본문보다 우선한다. 본문에서 `⟦옛 판 — …⟧` 가 붙은 문장은 인용하지 말 것.\n>\n" + blocks + "\n\n---\n\n")
    dst = src[:-3] + "_" + V + ".md"
    wr(dst, head + s)
    return dst


made = []
made.append(newver(R + "09_array_vs_flank/README_ARRAY_VS_FLANK.md",
    "> **정정 1(A9 · 옮겨 적기 오류)** : 인계문의 \"7 array 중 6 개가 flank 보다 +0.25~+0.41 · 종 공통\" 은 틀렸다. 이 README 의 그 문장은 **금화조 8 array 중 chr9 두 개를 뺀 6 개**에 대한 것이다.\n"
    "> CpG 합산으로 다시 잰 값(`00_catalog/CMD01_ARRAY_MINUS_FLANK.tsv`, array×시료) : 금화조 chr9 제외 6 array **+0.225~+0.415** · chr9 두 array −0.023~+0.101 · "
    "유인원 5 종 +0.184~+0.262(array 1 개씩) · 애기장대 9 array −0.112~+0.310(중앙 +0.008) · 제브라 31 array×6 개체 −0.281~+0.143(중앙 +0.003) · 생쥐 −0.21~−0.45(`10_mouse_5S`). **종 공통 명제가 아니다.**\n"
    "> **정정 2(A9 · 정의)** : §1 표의 \"금화조 10\" 은 출처 불명 — 원표 `AVF.tsv` 는 **8 array × 2 라이브러리**. 이 표의 array·flank '중앙' 은 중앙값 차라 CpG 합산 재실측과 ±0.02 다르다. 제브라 \"162\" 도 재실측 186(31×6).\n"
    "> **정정 3(A10 · 옛 해석)** : §3·§5 의 \"제브라는 발현과 맞는다(발현 좌위가 저메틸)\" 는 **철회**. chr18 이 꺼져 있을 알에서도 chr18 코딩부가 저메틸(−0.33~−0.53)이고, 저메틸은 **코딩부에만** 있다(`12_zebrafish_devseries/README_PARTS.md` · 인계문 §3-3). 성체 뇌 한 시점만 본 해석이었다.\n"
    "> **정정 4** : 애기장대 `CP116284.1` 은 **chr5**(`README_FIX_2026-09-13.md` 의 \"chr4 15.92 Mb\" 는 오기).",
    [("| 금화조 | 10 |", "정정 2"), ("## 3. 제브라피시 — 발현과 맞고", "정정 3"), ("발현 좌위(chr18)가 제 주변보다 0.09~0.13 낮다", "정정 3"),
     ("발현이 알려진 유일한 종에서, 발현 좌위가 주변 대비 저메틸이다", "정정 3"), ("특별한 것은 나머지 6 array 가 주변보다 +0.25~0.41 높다", "정정 1 : 금화조 8 중 chr9 제외 6"),
     ("제브라피시는 맞고(발현 좌위가 저메틸)", "정정 3")]))
made.append(newver(R + "00_catalog/TIER_STATUS.md",
    "> **정정 1(A10 · 해석 충돌)** : §5 \"세 종에서 저메틸 5S array 가 따로 나온다(수렴)\" 에서 **금화조 chr9 는 뺀다** — chr9 는 주변(flank)도 저메틸이라 array−flank 가 −0.023~+0.101(재실측)로 array 고유 저메틸이 아니다(`09_array_vs_flank`). "
    "남는 것은 생쥐(array−flank −0.21~−0.45)뿐이고 붉은털 FIXED 는 flank 대비를 재지 않아 **보류**. 따라서 '수렴' 주장은 쓰지 않는다.\n"
    "> **정정 2(A10)** : \"뇌 조직이라 다른 종(세포주)과 직접 비교하지 말 것\" — 제브라 ONT(PRJEB53713)도 **성체 뇌**다. 세포주는 유인원(GenomeArk) 쪽이다.\n"
    "> **정정 3(A9 · 정의)** : 붉은털 array 수 = TIER·지도 4 · 카탈로그·정본 5(NON_5S 1 포함) · 참조유닛 표 4 행(3 + NON_5S, NO_UNIT `CM111659.2_103988607` 없음). "
    "생쥐 `reads_placed 0` 은 'G8 read×유닛 표를 지난 read 수' 정의상 0(영역 추출 read 127 · 깊이 9.2 · junction 30). 정의 전부 `CMD01_CITE_DEFAULTS.tsv`.",
    [("## 5. 뜻밖의 수렴", "정정 1"), ("| 금화조 | chr9 두 array (372 bp 유닛)", "정정 1 : chr9 는 국소 환경"), ("다른 종(세포주)과 직접 비교하지 말 것", "정정 2")]))
made.append(newver(R + "09_noplacement/README_NOPLACEMENT.md",
    "> **정정 1(A10 · 옛 표현)** : \"배치 없이도 같은 답(값)이 나온다\" 는 과장이다. 재현되는 것은 **대비**(저메틸/고메틸 array 차 0.609–0.681)이고, "
    "**값 자체**는 같은 array 에서 경로 사이 0.000–0.068(중앙 0.023) 퍼진다. 동물 CpG 전용 모델에서 배치 경로와 −0.013~−0.019 일치는 유효. placement-free 이지 alignment-free 가 아니다.",
    [("**배치 없이도 같은 답이 나온다.**", "정정 1 : 같은 것은 대비"), ("## 3. 결과 — 금화조는 배치 없이도 같은 답이 나온다", "정정 1")]))
made.append(newver(R + "10_mouse_5S/README_MOUSE.md",
    "> **정정 1(A10)** : \"생쥐만 뇌 조직이고 나머지는 세포주/전체 개체\" — 제브라 ONT(PRJEB53713 대조 6 개체)도 **성체 뇌**다. 조직 차이로 교란되는 비교는 유인원(세포주)·금화조·애기장대 쪽이다.",
    [("생쥐만 뇌 조직이고 나머지는 세포주/전체 개체다", "정정 1 : 제브라도 성체 뇌")]))
made.append(newver(R + "12_zebrafish_devseries/README_DEVSERIES.md",
    "> **정정 1(A10 · 옛 값 잔존)** : Potok 알 chr18 은 **32M 판 −0.255 [−0.329, −0.181] r69**(옛 8M 판 −0.352), Potok 근육 chr18 은 **−0.098 [−0.148, −0.048] r145**(옛 −0.113). "
    "두 연구 알 구간(Jiang [−0.190, −0.034] · Potok [−0.329, −0.181])은 **가까스로 닿는다** — '구간이 겹친다' 는 옛 8M 기준 표현.\n"
    "> **정정 2(A10 · 옛 값 잔존)** : chr4 는 32M 판 Potok 알 **+0.038 [+0.032, +0.044]** · 근육 **+0.005 [+0.000, +0.010]**(옛 +0.037 [+0.025, +0.049] · +0.005 [−0.005, +0.015]).\n"
    "> **정정 3(A10)** : §5 \"난모세포 메틸 자료는 공개 자료에 없다\" 는 **틀렸다** — GSE175951 저자 가공 CpG 표(Stage III 난모세포 2 반복, `oocyte_GSE175951/README_OOCYTE.md`)와 GSE122695 생식세포 PBAT(`gsc_GSE122695/README_GSC.md`)가 있다. **Stage I–II** 난모세포 메틸만 없다.",
    [("**Potok −0.352**", "정정 1"), ("Potok  −0.352   sphere", "정정 1"), ("Potok 성체 근육은 +0.005 [−0.005, +0.015]", "정정 2"),
     ("반대로 Potok 알이 +0.037 [+0.025, +0.049]", "정정 2"), ("## 5. 난모세포 — chr4 에 필요한 시료는 공개 자료에 없다", "정정 3"),
     ("**메틸 자료는 없다.**", "정정 3"), ("chr18 부호가 같고(−0.112 vs −0.352) 구간이 겹친다", "정정 1")]))
made.append(newver(R + "12_zebrafish_devseries/README_SHORTREAD.md",
    "> **정정 1(A10)** : §2-2 끝 \"공개 자료에 쓸 수 있는 것이 없다\" 는 **틀렸다** — `PRJNA734348` 원시 read 는 우리 파이프라인에서 1–2% 만 붙지만 **저자 가공 CpG 표(GSE175951)는 정상**이라 사용한다. GSE122695(생식세포 PBAT)도 있다.\n"
    "> **정정 2** : 심장 chr4 는 read 단위 **+0.044**(옛 CpG 단위 +0.042).",
    [("공개 자료에 쓸 수 있는 것이 없다", "정정 1"), ("심장 +0.042", "정정 2")]))
for m in made:
    print("④ 새 판", m.replace(R, "data_v2/"))
