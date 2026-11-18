# -*- coding: utf-8 -*-
# Differential methylation at the array level: the tables and the tests, made only within one individual and one library.
#   The sampling unit is the read, and a read's value is the CpG methylation fraction it gives over that array, or locus.
#   1) The locus table: the whole-array, coding and spacer fractions per locus, with a read-level bootstrap 95% interval.
#   2) The tests: over the read values with at least MINC CpG sites, the fraction of variance the locus explains as eta squared, with a p value from 2,000 permutations of the locus label,
#            plus the key contrasts, zebra finch chr9 against the rest and zebrafish chr18 against chr4, with their differences and bootstrap 95% intervals.
#   Zebra finch chr2: the maternal and paternal arrays are byte-identical, so a haplotype cannot be assigned and the two are merged into one locus.
#   Zebrafish ONT: many chr4 arrays carry few reads, so only arrays with at least MINR reads count as a locus, against the single chr18 array.
# MIT License. Copyright (c) 2026 Hyungmin Kim
import os
import collections
import csv
import io
import numpy as np

rng = np.random.default_rng(0)
D = os.environ["RDNA_OTHERSP_DIR"] + "/data_v2/13_array_level_diff/"
T = "\t"
MINC = 10      # to be used as a read value, the read must carry at least 10 CpG sites at that locus
MINR = 20      # to enter the test as a locus, at least 20 reads are required
NB = 2000
NP = 2000

FINCH = {"CM109790.1_96427419_96511119": ("chr2A(mat/pat 배정불가)", "chr2", "배정불가"),
         "CM109751.1_96489183_96572883": ("chr2A(mat/pat 배정불가)", "chr2", "배정불가"),
         "CM109790.1_97181170_97425313": ("chr2B(mat/pat 배정불가)", "chr2", "배정불가"),
         "CM109751.1_97242934_97487077": ("chr2B(mat/pat 배정불가)", "chr2", "배정불가"),
         "CM109792.1_6514976_6563661": ("chr4_mat", "chr4", "mat"),
         "CM109753.1_6524260_6552707": ("chr4_pat", "chr4", "pat"),
         "CM109798.1_6937339_6952337": ("chr9_mat", "chr9", "mat"),
         "CM109759.1_7075841_7101562": ("chr9_pat", "chr9", "pat")}


def locus_of(sp, aid):
    core = aid.split("__")[1]
    if sp == "zfinch":
        return FINCH[core][0], FINCH[core][1]
    ch = "chr18" if core.startswith("CM111196.1") else "chr4"
    return ("chr18" if ch == "chr18" else "chr4:" + core.split("_", 1)[1]), ch


def load(sp):
    fn = D + "READ_ARRAY_PART_%s.tsv" % sp
    agg = collections.defaultdict(lambda: {"코딩부": [0, 0], "스페이서": [0, 0]})
    for r in csv.DictReader(io.open(fn, encoding="utf8"), delimiter=T):
        loc, ch = locus_of(sp, r["array_id"])
        k = (r["sample"], loc, ch, r["read_id"])
        x = agg[k][r["part"]]
        x[0] += int(r["n_calls"])
        x[1] += int(r["n_meth"])
    return agg


def frac(a):
    return a[:, 1].sum() / a[:, 0].sum() if a[:, 0].sum() else float("nan")


def boot(a):
    n = len(a)
    ix = rng.integers(0, n, size=(NB, n))
    return a[ix, 1].sum(1) / np.maximum(a[ix, 0].sum(1), 1)


tab = io.open(D + "ARRAY_LEVEL_TABLE.tsv", "w", encoding="utf8", newline="\n")
tab.write(T.join(["species", "sample", "locus", "chrom", "n_read", "whole", "whole_lo", "whole_hi", "coding", "coding_lo", "coding_hi",
                  "spacer", "spacer_lo", "spacer_hi", "read_median_whole"]) + "\n")
tst = io.open(D + "ARRAY_DIFF_TESTS.tsv", "w", encoding="utf8", newline="\n")
tst.write(T.join(["species", "sample", "test", "n_loci", "n_reads", "eta2", "perm_p", "contrast", "diff", "diff_lo", "diff_hi", "cohen_d"]) + "\n")

for sp in ("zfinch", "zebrafish"):
    agg = load(sp)
    bys = collections.defaultdict(lambda: collections.defaultdict(list))
    for (s, loc, ch, rid), v in agg.items():
        bys[s][(loc, ch)].append((v["코딩부"][0] + v["스페이서"][0], v["코딩부"][1] + v["스페이서"][1],
                                  v["코딩부"][0], v["코딩부"][1], v["스페이서"][0], v["스페이서"][1]))
    for s in sorted(bys):
        loci = bys[s]
        keep = {}
        for (loc, ch), rows in sorted(loci.items()):
            a = np.array(rows, dtype=float)
            w, c, sp_ = a[:, 0:2], a[:, 2:4], a[:, 4:6]
            bw, bc, bsp = boot(w), boot(c[c[:, 0] > 0]) if (c[:, 0] > 0).sum() > 1 else None, boot(sp_[sp_[:, 0] > 0]) if (sp_[:, 0] > 0).sum() > 1 else None
            rv = w[w[:, 0] >= MINC]
            med = np.median(rv[:, 1] / rv[:, 0]) if len(rv) else float("nan")
            def ci(b):
                return ("%.3f" % np.percentile(b, 2.5), "%.3f" % np.percentile(b, 97.5)) if b is not None else ("-", "-")
            tab.write(T.join([sp, s, loc, ch, str(len(a)), "%.3f" % frac(w), *ci(bw),
                              "%.3f" % frac(c), *ci(bc), "%.3f" % frac(sp_), *ci(bsp), "%.3f" % med]) + "\n")
            if len(rv) >= MINR:
                keep[(loc, ch)] = rv[:, 1] / rv[:, 0]
        # 2) eta squared and the permutation, over the read values against the locus label
        def eta2(groups):
            allv = np.concatenate(groups)
            gm = allv.mean()
            ssb = sum(len(g) * (g.mean() - gm) ** 2 for g in groups)
            sst = ((allv - gm) ** 2).sum()
            return ssb / sst if sst else float("nan")

        def perm(groups):
            obs = eta2(groups)
            allv = np.concatenate(groups)
            sizes = np.cumsum([len(g) for g in groups])[:-1]
            ge = 0
            for _ in range(NP):
                p = rng.permutation(allv)
                if eta2(np.split(p, sizes)) >= obs:
                    ge += 1
            return obs, (ge + 1) / (NP + 1)

        def contrast(a, b_, name, test, nl):
            ba = rng.choice(a, size=(NB, len(a))).mean(1)
            bb = rng.choice(b_, size=(NB, len(b_))).mean(1)
            dd = ba - bb        # independent samples, with the difference taken before any alignment
            sd = np.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b_) - 1) * b_.var(ddof=1)) / (len(a) + len(b_) - 2))
            tst.write(T.join([sp, s, test, str(nl), str(len(a) + len(b_)), "", "", name, "%+.3f" % (a.mean() - b_.mean()),
                              "%+.3f" % np.percentile(dd, 2.5), "%+.3f" % np.percentile(dd, 97.5), "%.2f" % ((a.mean() - b_.mean()) / sd)]) + "\n")

        groups = [v for v in keep.values()]
        if len(groups) >= 2:
            e, p = perm(groups)
            tst.write(T.join([sp, s, "모든 좌위 η²(read 값)", str(len(groups)), str(sum(len(g) for g in groups)), "%.3f" % e, "%.4f" % p, "", "", "", "", ""]) + "\n")
        if sp == "zfinch":
            c9 = [v for (l, ch), v in keep.items() if ch == "chr9"]
            ot = [v for (l, ch), v in keep.items() if ch != "chr9"]
            if c9 and ot:
                contrast(np.concatenate(c9), np.concatenate(ot), "chr9(mat+pat) − 나머지", "대비", len(keep))
            nc9 = [v for (l, ch), v in keep.items() if ch != "chr9"]
            if len(nc9) >= 2:
                e, p = perm(nc9)
                tst.write(T.join([sp, s, "chr9 뺀 좌위끼리 η²", str(len(nc9)), str(sum(len(g) for g in nc9)), "%.3f" % e, "%.4f" % p, "", "", "", "", ""]) + "\n")
            for ch in ("chr4", "chr9"):
                mm = keep.get((ch + "_mat", ch))
                pp = keep.get((ch + "_pat", ch))
                if mm is not None and pp is not None:
                    contrast(mm, pp, ch + " mat − pat", "하플로타입 대비", 2)
        else:
            c18 = keep.get(("chr18", "chr18"))
            c4 = [v for (l, ch), v in keep.items() if ch == "chr4"]
            if c18 is not None and c4:
                contrast(c18, np.concatenate(c4), "chr18 − chr4(read≥20 array 전부)", "대비", len(keep))
            if len(c4) >= 2:
                e, p = perm(c4)
                tst.write(T.join([sp, s, "chr4 array 끼리 η²", str(len(c4)), str(sum(len(g) for g in c4)), "%.3f" % e, "%.4f" % p, "", "", "", "", ""]) + "\n")
tab.close()
tst.close()
print("done")
