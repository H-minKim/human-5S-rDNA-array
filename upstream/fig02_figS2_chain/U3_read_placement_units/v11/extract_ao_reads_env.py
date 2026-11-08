# -*- coding: utf-8 -*-
# extract the stored accepted ARRAY_ONLY reads unchanged from the fastq (run on the node), header (MM/ML), sequence and quality unchanged
#   per sample, only the accepted reads of the one source that v7 used (GIAB > consensus2 (AOJV8) > below60 (AOJV5) > AOJV3)
#   the fastq is the file that job actually used on this node (srcmap, same path order as the original code)
#   Output: ~/aoreads_out/<sample>.<job>.<source>.fastq.gz  plus list ~/aoreads_out/manifest_<node>.tsv
import io
import os
import gzip
import socket
import collections

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
host = socket.gethostname()
I = H(os.environ.get("REPRO_IN", "~/repro_in"))


def samples(fn):
    return set(L.split(TAB)[0] for L in io.open(os.path.join(I, fn), encoding="utf8") if L.strip())


AGS = samples("AOJV_GIAB.tsv")
ACS = samples("AOJV8.tsv")
AFS = samples("AOJV5.tsv")


def cfg_of(s):
    return "AOJVG" if s in AGS else "AOJV8" if s in ACS else "AOJV5" if s in AFS else "AOJV3"


FN = {"AOJV3": "AOJV3.tsv", "AOJV5": "AOJV5.tsv", "AOJV8": "AOJV8.tsv", "AOJVG": "AOJV_GIAB.tsv"}
SAVED = collections.defaultdict(set)
for cfg, fn in FN.items():
    for L in io.open(os.path.join(I, fn), encoding="utf8"):
        r = L.split(TAB)
        if len(r) >= 2:
            SAVED[(cfg, r[0])].add(r[1])
O = H(os.environ.get("AOREADS", "~/aoreads_out"))
os.makedirs(O, exist_ok=True)
man = io.open(os.path.join(O, "manifest_%s.tsv" % host), "w", encoding="utf8", newline=NL)
for L in io.open(H("~/srcmap_%s.tsv" % host), encoding="utf8"):
    h, job, smp, lab, fq = L.rstrip(NL).split(TAB)[:5]
    if job not in FN or job != cfg_of(smp):
        continue
    want = SAVED.get((job, smp), set())
    if not want or fq == "-":
        continue
    out = os.path.join(O, "%s.%s.%s.fastq.gz" % (smp, job, lab))
    got = set()
    n = 0
    with gzip.open(out, "wt", encoding="utf8") as g:
        for p in fq.replace("~", H("~")).split(";"):
            with gzip.open(p, "rt", encoding="utf8") as fh:
                rec = []
                for x in fh:
                    rec.append(x)
                    if len(rec) == 4:
                        rid = rec[0][1:].split()[0]
                        if rid[:16] in want and rid[:16] not in got:
                            g.write("".join(rec))
                            got.add(rid[:16])
                            n += 1
                        rec = []
    man.write(TAB.join([host, job, smp, lab, fq, str(len(want)), str(n), os.path.basename(out)]) + NL)
man.close()
print(host, "done")
