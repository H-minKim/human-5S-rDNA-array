# -*- coding: utf-8 -*-
# fill the base qualities (QUAL) from fastq into the placement read SAM (run on the node)
#   the original ARRAY_ONLY code aligned from FASTA, so QUAL comes out as '*'. The alignment is kept and only QUAL is taken from the original read.
#   reverse (flag 16) records have SEQ reverse complemented, so QUAL is reversed too.
#   Check: fill only when the original read sequence (reverse complemented if reverse) and the SAM SEQ are identical character by character. Stop if they differ.
#   Arguments: input SAM, fastq paths (several joined with a semicolon), output SAM
import io
import sys
import gzip

NL = chr(10)
TAB = chr(9)
C = str.maketrans("ACGTNacgtn", "TGCANtgcan")
sam, fqs, out = sys.argv[1], sys.argv[2], sys.argv[3]
need = set()
for L in io.open(sam, encoding="utf8"):
    if not L.startswith("@"):
        need.add(L.split(TAB, 1)[0])
Q = {}
for fq in fqs.split(";"):
    with gzip.open(fq, "rt", encoding="utf8") as fh:
        rec = []
        for L in fh:
            rec.append(L.rstrip(NL))
            if len(rec) == 4:
                rid = rec[0][1:].split()[0]
                if rid in need and rid not in Q:
                    Q[rid] = (rec[1], rec[3])
                rec = []
n = filled = 0
o = io.open(out, "w", encoding="utf8", newline=NL)
for L in io.open(sam, encoding="utf8"):
    if L.startswith("@"):
        o.write(L)
        continue
    f = L.rstrip(NL).split(TAB)
    n += 1
    if f[10] == "*" and f[0] in Q:
        sq, qu = Q[f[0]]
        if int(f[1]) & 16:
            sq = sq.translate(C)[::-1]
            qu = qu[::-1]
        if sq.upper() != f[9].upper():
            sys.stderr.write("서열 불일치 %s - 중단%s" % (f[0], NL))
            raise SystemExit(1)
        f[10] = qu
        filled += 1
    o.write(TAB.join(f) + NL)
o.close()
print("%s 기록 %d · QUAL 채움 %d · fastq 에서 찾은 read %d/%d" % (sam.split("/")[-1], n, filled, len(Q), len(need)))
