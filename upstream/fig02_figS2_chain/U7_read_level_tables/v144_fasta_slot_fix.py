# -*- coding: utf-8 -*-
# v144 (2026-09-26) bug fix: deciding the slot from the file names in flank_final2.tar gets 53.2% wrong
#   finding: in 299 of 562 cases the slot in the file name (e.g. ..._1) and the slot of FLANK_OFFSET point to different haplotypes
#   evidence: the length of that fasta differs from the ref_len implied by its name, and matches the ref_len of the other slot of the same individual
#   affected: v143 (JUNCTION_SEQ_ANNOT, CpG density and Alu), qa_batch_B (CpG island statistics, 16-mer sharing), all on the sequence side
#   not affected: methylation outputs such as v142 (JUNCTION_FINE_METH) (they use only SAM+FLANK_OFFSET, not the fasta)
#   fix: link file to reference by sequence length, not by name. If lengths tie, resolve with the header flank5/flank3.
#   output FASTA_REF_MAP.tsv in the table folder (02_*) (canonical correspondence fasta name <-> ref_name), V144_MAP.txt in the log folder (04_*)
import io, os, tarfile
import pandas as pd
OUT = __import__("os").environ["RDNA_DATA_ROOT"] + "/v11_readvar"
TAR = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/12_flank/final/flank_final2.tar"
OFF = __import__("os").environ["RDNA_DATA_ROOT"] + "/for_paper/08_versions/v11/13_junction50k/mapping_table/FLANK_OFFSET.tsv"
LOG = io.open(OUT + "/04_logs/V144_MAP.txt", "w", encoding="utf8", newline="\n")
def P(s=""): print(s, flush=True); LOG.write(str(s) + "\n"); LOG.flush()

O = pd.read_csv(OFF, sep="\t", dtype={"slot": str})
by_sample = {}
for r in O.itertuples(index=False):
    by_sample.setdefault(r.sample, []).append(r)
def hd(txt):
    h = txt.split("\n", 1)[0]; m = {}
    for kv in h[1:].split("|"):
        if "=" in kv: k, v = kv.split("=", 1); m[k] = v
    return m, "".join(txt.split("\n")[1:])

rows = []; n_name_ok = n_fix = n_tie = n_none = 0
with tarfile.open(TAR) as tf:
    for m in tf.getmembers():
        if not m.name.endswith(".fa"): continue
        txt = tf.extractfile(m).read().decode()
        hdr, seq = hd(txt)
        base = os.path.basename(m.name)[:-3]
        tail = base.split("__")[-1]
        smp, slot = tail.rsplit("_", 1)
        smp = smp.replace("hiarr_", "").replace("arr_", "")
        cands = by_sample.get(smp, [])
        if not cands:
            rows.append((base, smp, slot, "", len(seq), "개체 없음")); n_none += 1; continue
        L = len(seq)
        hit = [c for c in cands if int(c.ref_len) == L]
        method = ""
        if len(hit) == 1:
            pick = hit[0]; method = "길이 일치"
        elif len(hit) > 1:
            # two haplotypes with the same length, resolved with the header flank5/flank3
            l5, l3 = int(hdr.get("flank5", -1)), int(hdr.get("flank3", -1))
            h2 = [c for c in hit if int(c.flank5_distal) == l5 or int(c.flank3_proximal) == l3]
            pick = h2[0] if h2 else hit[0]
            method = "길이 동률 → 헤더로 가림" if h2 else "길이 동률 → 첫 후보(주의)"
            n_tie += 1
        else:
            rows.append((base, smp, slot, "", L, "길이 맞는 참조 없음")); n_none += 1; continue
        named = [c for c in cands if str(c.slot) == slot]
        same_as_name = bool(named) and named[0].ref_name == pick.ref_name
        if same_as_name: n_name_ok += 1
        else: n_fix += 1
        rows.append((base, smp, slot, pick.ref_name, L, method + ("" if same_as_name else " · 이름과 다름")))
MAP = pd.DataFrame(rows, columns=["fasta", "sample", "slot_in_name", "ref_name", "fasta_len", "method"])
MAP.to_csv(OUT + "/02_tables/FASTA_REF_MAP.tsv", sep="\t", index=False)
P("[결과] fasta %d" % len(MAP))
P(" matched the name %d (%.1f%%), differed from the name and corrected %d (%.1f%%), length tied so the header was read %d, not found %d"
  % (n_name_ok, 100 * n_name_ok / len(MAP), n_fix, 100 * n_fix / len(MAP), n_tie, n_none))
P()
P("[확인] 고친 대응에서 ref_len 이 모두 맞는가")
chk = MAP[MAP.ref_name != ""].merge(O[["ref_name", "ref_len", "rcf", "ref_front", "canon_len"]], on="ref_name")
P("   대응된 %d · ref_len == fasta_len 인 비율 %.4f" % (len(chk), (chk.ref_len == chk.fasta_len).mean()))
P("   (참고) 이름 기반이었다면 맞았을 비율 %.4f" % (n_name_ok / len(MAP)))
LOG.close()
