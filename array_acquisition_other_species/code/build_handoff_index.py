# -*- coding: utf-8 -*-
# Builds the ledger of the newest files and the code index by measurement.
#   - LATEST_FILES.tsv : 최신판 항목마다 경로(글롭) · 파일 수 · 크기 · md5(여러 파일이면 md5 목록의 md5) · 수정일 · 설명 · 만든 코드
#   - LATEST_FILES_MANIFEST.tsv : 위 항목을 이루는 파일 하나하나의 md5
#   - CODE_INDEX.tsv : code/ 전 파일의 설명줄 · 수정일 · 상태(현행/옛판/판정기록/참고)
#   없는 경로가 있으면 멈춘다(조용한 누락 금지).
# Add an entry to SPEC when new work appears and run it again.
# MIT License · Copyright (c) 2026 Hungmin Kim
import glob
import hashlib
import io
import os
import sys
import time

R = "$RDNA_OTHERSP_DIR/"
OUT = R + "HANDOFF/"
TAB = "\t"
NL = "\n"

# (층, 경로 글롭(R 기준), 설명, 만든 코드)
SPEC = [
    # 00 카탈로그 · 판정표
    ("00_카탈로그", "data_v2/00_catalog/5S_array_catalog.tsv", "array 목록 170 행(조립본 16)", "array_catalog.py · cms_all.sh"),
    ("00_카탈로그", "data_v2/00_catalog/5S_array_structure.tsv", "array 구조표(처리한 12 조립본 · 77 array)", "array_struct.py"),
    ("00_카탈로그", "data_v2/00_catalog/5S_unit_status.tsv", "array 별 유닛 상태(09-13 수정 반영)", "make_status.py"),
    ("00_카탈로그", "data_v2/00_catalog/REFUNIT_PARAMS.tsv", "array 별 참조유닛 매개변수(주기·LO/HI·PAD·CTX) — 기본", "refunit_build.py · g8fix.py"),
    ("00_카탈로그", "data_v2/00_catalog/REFUNIT_PARAMS_MERGED.tsv", "위와 같되 제브라 56528778 을 477 bp(.k2) 구역값으로 — k2 지도용", "zf_all_arrays.py"),
    ("00_카탈로그", "data_v2/00_catalog/UNIT_REF_TABLE.tsv", "array 당 참조유닛 1 개(그룹핑 안 씀) 최종표 61 행", "unit_ref_table.py"),
    ("00_카탈로그", "data_v2/00_catalog/UNIT_TYPE_FINAL2.tsv", "유닛 종류 2 단계 관문 최종 판정(전 array 한 종류)", "unit_group_verify.py"),
    ("00_카탈로그", "data_v2/00_catalog/TIER_STATUS_FULL.tsv", "종별 Tier + 배치 건전성(깊이·관통·junction·mapq0·일치도)", "placement_quality.py"),
    ("00_카탈로그", "data_v2/00_catalog/5S_methylation_cross_species_summary.tsv", "종별 메틸 한 줄 요약(생쥐 행 없음 — 10_mouse_5S 참조)", "placement_meth_report.py"),
    ("00_카탈로그", "data_v2/00_catalog/PLACEMENT_METH_REPORT.tsv", "배치·메틸 실측 보고표", "placement_meth_report.py"),
    ("00_카탈로그", "data_v2/00_catalog/ZEBRAFISH_ALL_ARRAYS.tsv", "제브라 31 array 메틸·flank(ONT 성체 뇌)", "zf_all_arrays.py"),
    ("00_카탈로그", "data_v2/00_catalog/zebrafish_gene_types/zfgt.gene_hits.tsv", "제브라 cmsearch 5S 유전자 8,638 히트(유형 배정)", "zf_gene_types.py(노드)"),
    ("00_카탈로그", "data_v2/00_catalog/zf_regions.tsv", "제브라 array·flank 구역(GRCz12tu)", "수기(array_vs_flank 기준)"),
    ("00_카탈로그", "data_v2/00_catalog/dataset_access_status.tsv", "자료 접근 가능 여부(옥수수 막힘 등)", "수기"),
    # 01·02 정본 · 참조유닛 · 지도
    ("01_정본", "data_v2/01_canonical/*/*.fa", "array±20 kb 정본(flank-array-flank) + .unit.fa 참조유닛", "ape_arrays.sh · refunit_build.py · g8fix.py"),
    ("02_참조유닛", "data_v2/02_reference_units/*", "참조유닛 78 개 통합(09-13 수정 : 붉은털 7,523 · 애기장대 749 · 제브라 k2 477)", "refunit_build.py · g8fix.py"),
    ("02_지도", "data_v2/02_consensus_map/*/*", "G8 consensus map(참조 1 개짜리) — 유닛 좌표의 기준", "species_g8map.py · g8fix.py"),
    ("02_지도", "data_v2/02_consensus_map/G8MAP_SUM.tsv", "G8 지도 요약", "species_g8map.py"),
    # 03·04·05 read · BAM · bedMethyl
    ("04_배치BAM", "data_v2/04_bam/*/*.5S.bam", "배치 BAM 15(전 게놈 정렬 뒤 array 영역 추출 · MM/ML 보존)", "ape_align.sh · finch_select.sh · ath_finalize.sh · zf_extract.sh"),
    ("05_bedMethyl", "data_v2/05_bedmethyl/*/*", "modkit 자리별 메틸", "modkit(ape_to_d.sh · finch_finalize.sh · ath_finalize.sh)"),
    ("03_read×유닛메틸", "data_v2/03_reads/*/*", "read×유닛 메틸 표(G8 틀) + 유닛·array 요약", "units_matrix_g8.py"),
    ("06_유닛", "data_v2/06_units_from_canonical/*", "정본→유닛 분할(array_pos)", "canon_units.py"),
    ("06_유닛", "data_v2/06_units_from_reads/*", "read→유닛 분할(read_unit_seq, 변이 포함)", "read_units.py"),
    ("07_합본", "data_v2/07_consensus_from_reads/*", "read 합의 변이(깊이 35% · 엄격판 varlist_strict)", "unit_consensus.py · consensus_summary.py"),
    # 09 검정 · 종별
    ("09_검정", "data_v2/09_array_vs_flank/AVF.tsv", "array vs flank 100 kb(ONT) 원표", "array_vs_flank.py · region_meth.py"),
    ("09_검정", "data_v2/09_array_vs_flank/README_ARRAY_VS_FLANK_2026-09-21.md", "새 판 — 옛 README 는 HISTORY", "—"),
    ("09_검정", "data_v2/09_noplacement/*.tsv", "배치 없이 메틸(MM/ML 창 판별)", "meth_noplace.py"),
    ("09_검정", "data_v2/09_noplacement/README_NOPLACEMENT_2026-09-21.md", "새 판", "—"),
    ("09_검정", "data_v2/09_placement_quality/*", "배치 건전성 · 금화조 chr2 하플로타입 중복", "placement_quality.py · hap_assign.py"),
    ("09_검정", "data_v2/06_aux/lit_reference/*", "주기 직접측정·문헌대조·09-13 수정 기록(+fig3)", "ladder_all.py · period_direct2.py · lit_cmp.py · fig3.py"),
    ("09_검정", "data_v2/06_aux/ape_basecaller.tsv", "유인원 호출기 판본 표(판본 비교 금지 근거)", "ape_screen.sh"),
    ("10_생쥐", "data_v2/10_mouse_5S/*.ba*", "생쥐 영역 추출 BAM", "samtools view(영역)"),
    ("10_생쥐", "data_v2/10_mouse_5S/*.tsv", "생쥐 5S 메틸(무배치 창·read)", "meth_noplace.py · region_meth.py"),
    ("10_생쥐", "data_v2/10_mouse_5S/README_MOUSE_2026-09-21.md", "새 판", "—"),
    ("11_붉은털", "data_v2/11_rhesus_tier2/*", "붉은털원숭이 Tier2 자리별 메틸(bigWig)", "수기 추출"),
    # 12 제브라 short-read · 난모세포
    ("12_제브라", "data_v2/12_zebrafish_devseries/README_DEVSERIES_2026-09-21.md", "발생계열 새 판", "—"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/README_SHORTREAD_2026-09-21.md", "short-read 새 판", "—"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/README_PARTS.md", "코딩/스페이서 분해 문서", "—"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/README_ADULT_FLANK_OOGENESIS.md", "성체·flank·난자형성기 문서", "—"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/ADULT_STATS.tsv", "시료×좌위×부류(코딩/스페이서/flankL/flankR) read 단위 — 최신 종합", "gsc_perread_lr.py · adult_stats.py"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/ADULT_DIFFS_2026-09-21.tsv", "코딩·스페이서 − flank 차이 + read 수 · 인용 가능 표시(read<10 금지)", "adult_stats.py · cmd01_apply.py"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/ZFREAD_DEEP*.tsv", "array 전체 read 단위(깊이 32M 반영)", "bismark_read_meth.py · zf_deepen.sh"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/ZF_SHORTREAD_RUNS.tsv", "short-read 3 연구 전 run 목록 347", "ENA 조회"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/ENA_oocyte_egg_runs.tsv", "난모세포·알 관련 전 run 841", "ENA 조회"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/zf_regions_parts.tsv", "코딩/스페이서/flank 구역 조각 6,165", "build_zf_parts.py"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/gsc_GSE122695/*", "생식세포 저심도 PBAT(GSE122695) 군별 결과", "gsc_trim.py · gsc_align.sh · gsc_group.sh · gsc_perread.py · gsc_stats.py"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/oocyte_GSE175951/*", "Stage III 난모세포 가공 CpG(GSE175951) 결과", "find5S_danrer7.py · stem_5S.py"),
    ("12_제브라", "data_v2/12_zebrafish_devseries/oogenesis_5S_rna/*", "5S 유형 서열 · Locati 원문 · RNA 파일럿(메틸 아님)", "count5S.py"),
    # 문서
    ("문서", "data_v2/00_catalog/TIER_STATUS_2026-09-21.md", "Tier 문서 새 판", "—"),
    ("00_카탈로그", "data_v2/00_catalog/CMD01_*.tsv", "CMD-01 실측표 : array−flank · array 수 정의 · 인용 기본값 · 극성 · 보노보 판본별", "cmd01_measure.py · cmd01_apply.py · cmd01_bonobo_basecaller.py"),
    ("문서", "data_v2/00_catalog/UNIT_TYPES_README.md", "유닛 종류 판정 문서(결론: 그룹핑 안 씀)", "—"),
    ("문서", "data_v2/00_catalog/PLACEMENT_METH_README.md", "배치·메틸 방법 문서", "—"),
    ("문서", "data_v2/00_catalog/ZEBRAFISH_ALL_ARRAYS.md", "제브라 31 array 문서(§4 발현 해석은 옛판 — HISTORY §2)", "—"),
    ("문서", "data_v2/06_units_from_canonical/README_UNITS.md", "유닛 분할 문서", "—"),
    ("문서", "data_v2/07_consensus_from_reads/README_CONSENSUS.md", "07 합본 문서", "—"),
    # array 단위 차등 메틸(주 주장)
    ("13_array단위", "data_v2/13_array_level_diff/README_ARRAY_LEVEL_DIFF.md", "주 주장 문서 : array 단위 차등 메틸(금화조·제브라)", "—"),
    ("13_array단위", "data_v2/13_array_level_diff/ARRAY_*.tsv", "좌위별 전체·코딩·스페이서 + 부트스트랩 · η²·순열·대비", "cmd02_stats.py"),
    ("13_array단위", "data_v2/13_array_level_diff/ZF_CHR4_ARRAY_REPRO*.tsv", "제브라 chr4 array × 개체 · 개체 간 순위 상관", "cmd02_reprod.py"),
    ("13_array단위", "data_v2/13_array_level_diff/WGBS_*.tsv", "제브라 short-read 시료별 chr18−chr4 + read 단위 입력", "cmd02_wgbs.py"),
    ("13_array단위", "data_v2/13_array_level_diff/REFUNIT_IDENTITY.tsv", "참조유닛 스페이서·유전자 일치도(원인 후보 보조)", "cmd02_seqident.py"),
    ("13_array단위", "data_v2/13_array_level_diff/READ_ARRAY_PART_*.tsv", "read×array×코딩/스페이서 CpG 호출(ONT MM/ML)", "cmd02_extract.py"),
    ("13_array단위", "data_v2/13_array_level_diff/REGIONS_*.tsv", "array 구간·코딩부 구간", "cmd02_regions.py"),
    # 그림 폴더(2026-10-01 · 판3)
    ("그림", "FIGURES/README_figure_rules.md", "그림 규칙 : 세션 전용 폴더 저장 · 이름은 세션 그대로 · 태그(타종/othersp)는 보낼 때만", "—"),
    ("결과대장", "RESULTS/README_folder_rules.md", "결과 대장 판 규칙(옛 판은 _previous/)", "—"),
    ("결과대장", "RESULTS/latest_2026-10-29/*", "결과 28 행 × 산출·코드·입력·재현 + read 수준 접근번호(58 flowcell→52 run)", "run_repro.py · make_ledger.py · fetch_acc.py"),
    # 노드 정리 백업(2026-09-27)
    ("노드정리", "data_v2/00_catalog/NODE_CLEANUP_INVENTORY_2026-09-27.md", "노드 재다운로드 목록 : URL·계정·실측 속도·용도", "—"),
    ("노드정리", "data_v2/00_catalog/NODE_REBUILD_RECIPE_2026-09-27.md", "되살리는 법 : 금화조 조립본 레시피(md5 증명)·회수분·남은 단계", "—"),
    ("노드정리", "data_v2/00_catalog/NODE_MANIFEST_2026-09-27/*.tsv", "노드 3 대 파일 목록(크기·수정일·경로 · 24,504 파일)", "—"),
    ("노드정리", "data_v2/00_catalog/NODE_SCRIPTS_2026-09-27/*/scripts.tgz", "노드 실행 스크립트 회수본(113 파일 · 116 KB)", "—"),
    ("06_보조", "data_v2/06_aux/bonobo_read_basecaller.tsv", "보노보 read→호출기 판본 대응표(노드에서 회수)", "cmd01_bonobo_basecaller.py"),
    # 유인원 유닛 변이표( 제공분)
    ("14_유인원유닛변이", "data_v2/14_ape_unit_variants/README_orientation_convention_and_tables.md", "방향 규약(array_pos 1 = 사람과 같은 쪽 · flank 상동 실측) + 표 설명", "cmd03_flank_ortho.py"),
    ("14_유인원유닛변이", "data_v2/14_ape_unit_variants/APE_*.tsv", "유인원 6 종 read×유닛 · 유닛자리 합의 변이(종 정본 대비) · 좌표", "cmd03_ape_unit_tables.py"),
    ("14_유인원유닛변이", "data_v2/14_ape_unit_variants/APE_REFERENCE_UNITS.fa", "종별 정본 유닛 FASTA(+사람 CHM13)", "cmd03_ape_unit_tables.py"),
    ("14_유인원유닛변이", "data_v2/14_ape_unit_variants/FLANK_ORTHOLOGY.tsv", "flank 21-mer 상동 대조(방향 판정 근거)", "cmd03_flank_ortho.py"),
]

# code/ 상태 — 기재 없으면 '현행'
CODE_STATUS = {
    "g8map_multiref.py": "철회(다중참조 지도 — 좌표 변화 0 · 쓰지 않음)",
    "multiref_map.py": "철회(다중참조)",
    "unit_group_call.py": "판정기록(유닛 그룹 판정 — 결론: 그룹핑 안 씀)",
    "unit_group_verify.py": "판정기록(2 단계 관문 · UNIT_TYPE_FINAL2 생성)",
    "unit_types.py": "판정기록", "unit_types_gate.py": "판정기록", "unit_types_merge.py": "판정기록",
    "group_holdout.py": "판정기록",
    "species_map.py": "옛판(G6 이전 지도 → species_g8map.py)",
    "species_map2.py": "옛판(G6 지도 → species_g8map.py)",
    "units_matrix.py": "옛판(G6 틀 → units_matrix_g8.py)",
    "zf_units.py": "옛판(초기 제브라 유닛)",
    "bismark_region_meth.py": "참고(CpG 단위 — 결론은 read 단위 스크립트로)",
    "zf_qc_runs.sh": "옛판(단가닥 변환검사 → zf_qc_strand.sh)",
    "zf_qc_eggs2.sh": "옛판(단가닥 변환검사 → zf_qc_strand.sh)",
    "zf_parts_sum.py": "옛판(ZFPARTS 정리 → adult_stats.py)",
    "gsc_perread.py": "현행(flank 합침판) — flank 좌우는 gsc_perread_lr.py",
    "count5S.py": "참고(RNA 5S 유형 계수 — 메틸 아님)",
    "cmd01_measure.py": "현행", "cmd01_apply.py": "현행",
    "cmd01_prompt.py": "현행", "cmd01_bonobo_basecaller.py": "현행",
    "build_handoff_index.py": "현행(인계 대장 생성기)",
    "cmd02_regions.py": "현행", "cmd02_extract.py": "현행",
    "cmd02_stats.py": "현행", "cmd02_reprod.py": "현행",
    "cmd02_wgbs.py": "현행", "cmd02_seqident.py": "현행",
    "cmd03_ape_unit_tables.py": "현행",
    "cmd03_flank_ortho.py": "현행",
}


def md5f(p):
    h = hashlib.md5()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def desc_of(p):
    try:
        lines = io.open(p, encoding="utf8", errors="replace").read().splitlines()[:12]
    except Exception:
        return ""
    for L in lines:
        s = L.strip()
        if s.startswith("#!") or "coding" in s[:20] or not s.startswith("#"):
            continue
        s = s.lstrip("#").strip()
        if s:
            return s[:140]
    return ""


os.makedirs(OUT, exist_ok=True)
missing = []
lf = io.open(OUT + "LATEST_FILES.tsv", "w", encoding="utf8", newline=NL)
mf = io.open(OUT + "LATEST_FILES_MANIFEST.tsv", "w", encoding="utf8", newline=NL)
lf.write(TAB.join(["layer", "path_glob", "n_files", "bytes", "md5_or_md5ofmd5", "latest_mtime", "description", "made_by"]) + NL)
mf.write(TAB.join(["layer", "path", "bytes", "md5", "mtime"]) + NL)
for layer, g, d, code in SPEC:
    fs = sorted(p for p in glob.glob(R + g) if os.path.isfile(p))
    if not fs:
        missing.append(g)
        continue
    sums = []
    tot = 0
    mt = 0
    for p in fs:
        m = md5f(p)
        s = os.path.getsize(p)
        t = os.path.getmtime(p)
        sums.append(m)
        tot += s
        mt = max(mt, t)
        mf.write(TAB.join([layer, p[len(R):], str(s), m, time.strftime("%Y-%m-%d %H:%M", time.localtime(t))]) + NL)
    agg = sums[0] if len(sums) == 1 else hashlib.md5("".join(sums).encode()).hexdigest()
    lf.write(TAB.join([layer, g, str(len(fs)), str(tot), agg, time.strftime("%Y-%m-%d %H:%M", time.localtime(mt)), d, code]) + NL)
lf.close()
mf.close()

cf = io.open(OUT + "CODE_INDEX.tsv", "w", encoding="utf8", newline=NL)
cf.write(TAB.join(["file", "mtime", "status", "description"]) + NL)
for p in sorted(glob.glob(R + "code/*")):
    b = os.path.basename(p)
    cf.write(TAB.join([b, time.strftime("%Y-%m-%d", time.localtime(os.path.getmtime(p))),
                       CODE_STATUS.get(b, "현행"), desc_of(p)]) + NL)
cf.close()
if missing:
    sys.stderr.write("없는 경로 : %r\n" % missing)
    sys.exit(1)
print("LATEST 항목 %d · CODE %d" % (len(SPEC), len(glob.glob(R + "code/*"))))
