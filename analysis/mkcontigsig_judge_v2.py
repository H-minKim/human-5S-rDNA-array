# Significance call v2: criterion 3 changes from "median over read pairs >= 3" to "fraction of disagreeing pairs (score < 0) <= 0.2". In a large contig a random pair of reads overlaps only a little, so a low median comes from too little overlap rather than from disagreement.
import json, sys, collections
sys.stdout.reconfigure(encoding="utf-8")
def judge(o):
    why = []
    if o["n_mol"] < 2: why.append("독립 read <2")
    if o["n_q20"] < 2: why.append("변이 염기가 Q>=20 인 read %d 개(변이 Q %s · 같은 read 같은 자리 Q 중앙과의 차 %s)" % (o["n_q20"], o["q_alt"], o["dq_med"]))
    if o["pair_neg"] is None or o["pair_neg"] > 0.2: why.append("read 끼리 어긋남(점수<0 쌍 비율 %s)" % o["pair_neg"])
    return ("유의" if not why else "유의 아님"), "; ".join(why)
if __name__ == "__main__":
    for tag in sys.argv[1:]:
        O = json.load(open("mk2/%s_contigsig.json" % tag, encoding="utf-8"))
        for o in O: o["sig2"], o["why2"] = judge(o)
        json.dump(O, open("mk2/%s_contigsig_v2.json" % tag, "w", encoding="utf-8"), ensure_ascii=False, indent=0)
        print("==", tag, "contig", len(O))
        c = collections.Counter((o["cls"], o["sig2"]) for o in O)
        for k in sorted(c): print("   %-34s %-6s %d" % (k + (c[k],)))
        for o in O:
            if o["cls"][0] != "A" or o["sig2"] != "유의":
                print("    [%s · %s] %-6s %s read %d %s · 변이 Q %s · Q20 read %d · 같은자리 Q차 중앙 %s · 어긋난 쌍 %s · %s" % (o["cls"][:6], o["sig2"], o["variant"], o["mzone"], o["n_read"], o["plat"], o["q_alt"][:10], o["n_q20"], o["dq_med"], o["pair_neg"], o["why2"] or "-"))
