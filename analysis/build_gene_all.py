# -*- coding: utf-8 -*-
"""전체 hits read(배치기 안 거침) 판 — 유전자부 β 층 라벨 + 풀링 유닛 β 이봉성 요약. 노드에서 실행."""
import io, os, collections
NL = chr(10); TB = chr(9)
D = os.path.expanduser("~/cells")
os.makedirs(os.path.join(D, "strata"), exist_ok=True)
for lab, ds in (("K562all", "K562all"), ("GM12878all", "GM12878all")):
    f = io.open(os.path.join(D, "UNIT2_%s.tsv" % lab), encoding="utf8")
    ix = dict((k, i) for i, k in enumerate(f.readline().rstrip(NL).split(TB)))
    g = collections.OrderedDict()
    for l in f:
        x = l.rstrip(NL).split(TB)
        if x[ix["beta_gene"]] != "" and int(x[ix["n_cpg_gene"]]) >= 5:
            g[x[ix["unit_key"]]] = "low" if float(x[ix["beta_gene"]]) < 0.5 else "high"
    with io.open(os.path.join(D, "strata", "%s_gene.tsv" % ds), "w", encoding="utf8", newline=NL) as o:
        o.write("unit_key" + TB + "label" + NL)
        for k, v in g.items():
            o.write(k + TB + v + NL)
    print(lab, "유전자부 층", dict(collections.Counter(g.values())))
