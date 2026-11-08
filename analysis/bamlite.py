# -*- coding: utf-8 -*-
"""순수 파이썬 BAM 읽기 (pysam 없음) — BGZF 는 연결된 gzip 멤버라 gzip 모듈로 풀린다.
필요한 것만: 이름·flag·ref·pos·CIGAR·SEQ·MM·ML.
"""
import gzip, struct

_SEQ = "=ACMGRSVTWYHKDBN"
_CIG = "MIDNSHP=X"


def _aux(b, i, end):
    tags = {}
    while i < end:
        tag = b[i:i + 2].decode(); t = chr(b[i + 2]); i += 3
        if t in "cC":
            v = b[i]; i += 1
        elif t in "sS":
            v = struct.unpack_from("<h" if t == "s" else "<H", b, i)[0]; i += 2
        elif t in "iI":
            v = struct.unpack_from("<i" if t == "i" else "<I", b, i)[0]; i += 4
        elif t == "f":
            v = struct.unpack_from("<f", b, i)[0]; i += 4
        elif t == "A":
            v = chr(b[i]); i += 1
        elif t in "ZH":
            j = b.index(0, i); v = b[i:j].decode(); i = j + 1
        elif t == "B":
            st = chr(b[i]); n = struct.unpack_from("<i", b, i + 1)[0]; i += 5
            fmt = {"c": "b", "C": "B", "s": "h", "S": "H", "i": "i", "I": "I", "f": "f"}[st]
            sz = struct.calcsize(fmt)
            v = list(struct.unpack_from("<%d%s" % (n, fmt), b, i)); i += n * sz
        else:
            raise ValueError("aux type " + t)
        tags[tag] = v
    return tags


def read_bam(path):
    with gzip.open(path, "rb") as f:
        data = f.read()
    assert data[:4] == b"BAM\x01"
    i = 4
    l_text = struct.unpack_from("<i", data, i)[0]; i += 4 + l_text
    n_ref = struct.unpack_from("<i", data, i)[0]; i += 4
    refs = []
    for _ in range(n_ref):
        ln = struct.unpack_from("<i", data, i)[0]; i += 4
        name = data[i:i + ln - 1].decode(); i += ln
        i += 4
        refs.append(name)
    while i < len(data):
        bs = struct.unpack_from("<i", data, i)[0]
        rec = data[i + 4:i + 4 + bs]; i += 4 + bs
        refid, pos, lrn, mapq, _bin, ncig, flag, lseq = struct.unpack_from("<iiBBHHHi", rec, 0)
        o = 32
        name = rec[o:o + lrn - 1].decode(); o += lrn
        cig = []
        for k in range(ncig):
            v = struct.unpack_from("<I", rec, o + 4 * k)[0]
            cig.append((_CIG[v & 15], v >> 4))
        o += 4 * ncig
        sb = rec[o:o + (lseq + 1) // 2]; o += (lseq + 1) // 2
        seq = "".join(_SEQ[(sb[k // 2] >> (4 if k % 2 == 0 else 0)) & 15] for k in range(lseq))
        o += lseq
        tags = _aux(rec, o, len(rec))
        yield dict(name=name, flag=flag, ref=refs[refid] if refid >= 0 else None,
                   pos=pos, cigar=cig, seq=seq, tags=tags)
