# -*- coding: utf-8 -*-
# per node and per sample, the fastq actually used for placement and its data source (run on the node)
#   chosen with the same path order as rjfreq2 / aoplace (the first file that exists)
#   source labels: HPRC_ONT (s5s_ont = UCSC HPRC submission, s5s_new = HPRC working new 44)
#               1kGP_ONT (kgp_ont/out) · GIAB_ONT (giab5s_seq · ONT Open Data)
#               GIAB_SOMATIC_ONT (hg008/out), PACBIO_HIFI (cancer5s/hits, HG008-T and COLO829 are in the HiFi manifest)
import io
import os
import sys
import glob
import socket

NL = chr(10)
TAB = chr(9)
H = os.path.expanduser
ORDER = [("~/s5s_ont/hits/%s.hits.fastq.gz", "HPRC_ONT_submissions"),
         ("~/s5s_new/hits/%s.hits.fastq.gz", "HPRC_ONT_working"),
         ("~/kgp_ont/out/%s.array.fastq.gz", "1kGP_ONT"),
         ("~/cancer5s/hits/%s.hits.fastq.gz", "PACBIO_HIFI"),
         ("~/hg008/out/%s.hits.fastq.gz", "GIAB_SOMATIC_ONT")]


def pick(smp):
    for pat, lab in ORDER:
        q = H(pat % smp)
        if os.path.exists(q) and os.path.getsize(q) > 1000:
            return q, lab
    ch = sorted(x for x in glob.glob(H("~/giab5s_seq/out/%s_*.hits.fastq.gz" % smp))
                + glob.glob(H("~/hg008/out/%s_*.hits.fastq.gz" % smp)) if os.path.getsize(x) > 1000)
    if ch:
        lab = "GIAB_ONT" if "giab5s_seq" in ch[0] else "GIAB_SOMATIC_ONT"
        return ";".join(ch), lab
    return "-", "NONE"


host = socket.gethostname()
o = io.open(H("~/srcmap_%s.tsv" % host), "w", encoding="utf8", newline=NL)
for job, lst in (("JX", "~/rp_j.txt"), ("JX_C2", "~/rp_j2.txt"), ("AOJV3", "~/rp_a3.txt"),
                 ("AOJV5", "~/rp_a5.txt"), ("AOJV8", "~/rp_c2here.txt"), ("AOJVG", "~/rp_g.txt")):
    p = H(lst)
    if not os.path.exists(p):
        continue
    for s in sorted(set(x.strip() for x in io.open(p, encoding="utf8") if x.strip())):
        fq, lab = pick(s)
        o.write(TAB.join([host, job, s, lab, fq.replace(H("~"), "~")]) + NL)
o.close()
print(host, "done")
