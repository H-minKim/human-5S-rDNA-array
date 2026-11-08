# -*- coding: utf-8 -*-
r"""The rulers, as a function instead of a paragraph.

Why this file exists
--------------------
`CONVENTIONS.md` has carried a table of rulers since the archive was started, as its very
That is the lesson this file is the answer to: a ruler that lives in prose is followed when
somebody remembers to look it up, and a ruler that lives in a function is followed by anyone
who calls the function. `CONVENTIONS.md` keeps the reasons; this keeps the rulers.

    from rulers import content_md5, same          # checks/ is on sys.path when run from here

Each ruler below states the false conclusion it prevents. That sentence is the point of the
ruler, not decoration: every one of them was written after a real wrong answer.
"""
import gzip
import hashlib
import io
import os
import re

__all__ = ['content_md5', 'same', 'RULERS', 'ruler_for']


class RulerError(Exception):
    """Raised instead of returning a sentinel.

    A checker that turns an unreadable file into the string `err` will later compare `err`
    with `err` and report a match. That happened here on 2026-10-22: two failures were
    counted as an agreement. Never return a value that can equal another failure.
    """


def _md5_bytes(b):
    return hashlib.md5(b).hexdigest()


#  ---------------------------------------------------------------- the rulers


def _ruler_gz(path):
    """`.gz`: hash the decompressed content, never the file.

    False conclusion prevented: "the table changed". gzip writes the current time into its
    header, so compressing identical bytes twice gives two different file md5s. Measured
    2026-10-01 on `MODEL_SITE_V11.tsv.gz` (0071d1a855e9 against 00111bc77a0f, same content)
    and again 2026-10-23 on three P10 tables of 12 to 14 million lines.
    """
    h = hashlib.md5()
    with gzip.open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


class _Digest(object):
    """A file-like sink that hashes what is written to it and keeps nothing."""

    def __init__(self):
        self._h = hashlib.md5()

    def write(self, text):
        self._h.update(text.encode('utf-8'))
        return len(text)

    def hexdigest(self):
        return self._h.hexdigest()


def _ruler_parquet(path):
    """`.parquet`: hash `to_csv(index=False, lineterminator='\n')`, never the file.

    False conclusion prevented: "the table changed". Row-group layout, compression and
    writer version all change the bytes while the table is the same.


    1. "the table changed between my machine and the node". `to_csv` defaults its line
       terminator to `os.linesep`, so the same frame hashes differently on Windows
       (`\r\n`) and on Linux (`\n`). A ruler that depends on the operating system
       cannot compare a local copy with a node copy, which is most of what this ruler
       is for. The terminator is therefore named here and must stay named.
    2. "the ruler failed, so something is wrong with the file". Building the whole CSV
       as one string raised MemoryError on a 60 MB parquet of 733,288 rows, and passed
       on a second run with more memory free. A measurement that depends on free memory
       is not a ruler. The frame is streamed into the digest instead, so the result does
       not depend on how much memory happens to be available.
    """
    try:
        import pandas as pd
    except ImportError as e:                                    # pragma: no cover
        raise RulerError('pandas is needed to read %s' % path) from e
    sink = _Digest()
    pd.read_parquet(path).to_csv(sink, index=False, lineterminator='\n',
                                 chunksize=50000)
    return sink.hexdigest()


def _ruler_markdown(path):
    """`.md`: hash the raw bytes, with only the line endings normalised.

    False conclusion prevented: "the two documents are the same" when their section headings
    differ. Markdown headings all start with `#`, so the text ruler below -- written to drop a
    `#` provenance header from a table -- silently removes every heading from a document. A
    minimal case: two files whose only difference is `# Methods` against `# Data availability`
    hashed identically under that ruler, and `same()` returned True. Line endings are still
    normalised, because this archive is written on Windows and a `CRLF` would otherwise move
    the digest without the content moving.
    """
    b = io.open(path, 'rb').read().replace(b'\r\n', b'\n')
    return hashlib.md5(b).hexdigest()


def _ruler_text(path):
    """`.tsv` / `.csv` / `.txt`: hash the body, with `#` header lines removed.

    False conclusion prevented: "the table changed" when only a provenance header was added
    or a timestamp line rewritten. Compare bodies; compare headers separately and on purpose.
    """
    h = hashlib.md5()
    with io.open(path, 'rb') as f:
        for line in f:
            if line.startswith(b'#'):
                continue
            h.update(line.replace(b'\r\n', b'\n'))
    return h.hexdigest()


_SVG_ID = re.compile(rb'id="([^"]+)"')
_SVG_DATE = re.compile(rb'<dc:date>[^<]*</dc:date>')
_SVG_ISO = re.compile(rb'\d{4}-\d{2}-\d{2}T[\d:.]+')


def _ruler_svg(path):
    """`.svg`: normalise the renderer's random ids and its date, then hash -- and only
    compare two files from the same renderer.

    False conclusion prevented: "the figure changed". matplotlib mints a fresh `clipPath`
    id on every render, so two renders of one figure never match byte for byte. Normalising
    is valid *within* one renderer; across renderers an svg says nothing, and the pixel
    ruler is the one to use.
    """
    b = io.open(path, 'rb').read()
    b = _SVG_DATE.sub(b'<dc:date></dc:date>', b)
    b = _SVG_ISO.sub(b'__t__', b)
    ids = sorted(set(_SVG_ID.findall(b)), key=lambda x: (-len(x), x))
    for i, old in enumerate(ids):
        b = b.replace(old, ('__id%d__' % i).encode())
    return '%s+ids%d' % (_md5_bytes(b), len(ids))


def _ruler_pixels(path):
    """`.png` / `.pdf`: hash the rendered pixels, never the file.

    False conclusion prevented: both directions. A PDF carries `/CreationDate`, so identical
    drawings differ; and a PNG written by a different writer differs byte for byte while the
    picture is the same. Pixels are the only ruler that answers "does it look the same".
    """
    try:
        import fitz
    except ImportError as e:                                    # pragma: no cover
        raise RulerError('PyMuPDF (fitz) is needed for the pixel ruler on %s' % path) from e
    doc = fitz.open(path)
    h = hashlib.md5()
    for page in doc:
        h.update(page.get_pixmap().samples)
    doc.close()
    return h.hexdigest()


def _ruler_bytes(path):
    """Anything else: the file's own md5.

    Use it only for files that no tool rewrites -- source code, a `.fa`, a `.bam` you did
    not make. For anything a program emits, find the ruler above that fits, or add one here
    and say in `CONVENTIONS.md` why.
    """
    h = hashlib.md5()
    with io.open(path, 'rb') as f:
        for chunk in iter(lambda: f.read(1 << 20), b''):
            h.update(chunk)
    return h.hexdigest()


RULERS = (
    (('.gz',), 'decompressed content', _ruler_gz),
    (('.parquet',), "to_csv(index=False, lineterminator='\\n')", _ruler_parquet),
    (('.md',), 'raw bytes, CRLF normalised', _ruler_markdown),
    (('.tsv', '.csv', '.txt'), 'body without # header lines', _ruler_text),
    (('.svg',), 'ids and date normalised (one renderer only)', _ruler_svg),
    (('.png', '.pdf'), 'rendered pixels', _ruler_pixels),
)


def ruler_for(path):
    """Name the ruler that applies, without reading the file.

    Print this beside every number you report. A comparison whose ruler is not stated cannot
    be lined up against anyone else's comparison, which is how two correct measurements come
    to look like a disagreement.
    """
    ext = os.path.splitext(path)[1].lower()
    for exts, name, _ in RULERS:
        if ext in exts:
            return name
    return 'file bytes'


def content_md5(path):
    """The md5 of `path` under the ruler that fits its kind. Raises `RulerError`, never a
    sentinel."""
    if not os.path.exists(path):
        raise RulerError('no such file: %s' % path)
    ext = os.path.splitext(path)[1].lower()
    for exts, _, fn in RULERS:
        if ext in exts:
            return fn(path)
    return _ruler_bytes(path)


def same(a, b):
    """True when `a` and `b` agree under their ruler. Both must be the same kind of file,
    and either one failing to read is an exception rather than a mismatch."""
    ea, eb = os.path.splitext(a)[1].lower(), os.path.splitext(b)[1].lower()
    if ea != eb:
        raise RulerError('different kinds of file, so no ruler covers both: %s vs %s' % (ea, eb))
    return content_md5(a) == content_md5(b)


#  ---------------------------------------------------------------- self-test


def _selftest():
    """Show that each ruler separates what it should and joins what it should.

    A ruler that said `True` for everything would pass every comparison in this archive and
    be worse than none. Each case below therefore has a positive and a negative half.
    """
    import shutil
    import tempfile
    import time
    d = tempfile.mkdtemp(prefix='rulers_')
    ok = []
    try:
        #  ⑴ .gz -- same content, two compressions a second apart: file md5 differs, ruler agrees
        p1, p2 = os.path.join(d, 'a.tsv.gz'), os.path.join(d, 'b.tsv.gz')
        payload = b'col\n' + b''.join(b'%d\n' % i for i in range(5000))
        with gzip.open(p1, 'wb') as f:
            f.write(payload)
        time.sleep(1.1)
        with gzip.open(p2, 'wb') as f:
            f.write(payload)
        ok.append(('gz: file bytes differ', _ruler_bytes(p1) != _ruler_bytes(p2)))
        ok.append(('gz: ruler agrees', same(p1, p2)))
        with gzip.open(p2, 'wb') as f:
            f.write(payload + b'9999\n')
        ok.append(('gz: ruler separates changed content', not same(p1, p2)))

        #  ⑵ text -- a # header added: body ruler agrees, file bytes do not
        t1, t2 = os.path.join(d, 'a.tsv'), os.path.join(d, 'b.tsv')
        io.open(t1, 'wb').write(b'x\ty\n1\t2\n')
        io.open(t2, 'wb').write(b'# built 2026-10-23\nx\ty\n1\t2\n')
        ok.append(('tsv: file bytes differ', _ruler_bytes(t1) != _ruler_bytes(t2)))
        ok.append(('tsv: ruler agrees', same(t1, t2)))
        io.open(t2, 'wb').write(b'# built\nx\ty\n1\t3\n')
        ok.append(('tsv: ruler separates changed body', not same(t1, t2)))

        #  ⑶ markdown -- headings begin with # and must survive the ruler. This is the case
        #     that the text ruler failed for three weeks: two documents differing only in
        #     their section headings hashed identically and same() returned True.
        m1, m2 = os.path.join(d, 'a.md'), os.path.join(d, 'b.md')
        io.open(m1, 'wb').write(b'# Methods\n\nThe cut is 0.508.\n')
        io.open(m2, 'wb').write(b'# Data availability\n\nThe cut is 0.508.\n')
        ok.append(('md: ruler separates changed heading', not same(m1, m2)))
        io.open(m2, 'wb').write(b'# Methods\r\n\r\nThe cut is 0.508.\r\n')
        ok.append(('md: ruler ignores CRLF', same(m1, m2)))
        io.open(m2, 'wb').write(b'# Methods\n\nThe cut is 0.509.\n')
        ok.append(('md: ruler separates changed body', not same(m1, m2)))

        #  ⑷ failures raise -- they must not compare equal to each other
        try:
            same(os.path.join(d, 'nope1.tsv'), os.path.join(d, 'nope2.tsv'))
            ok.append(('missing files raise', False))
        except RulerError:
            ok.append(('missing files raise', True))
        try:
            same(t1, p1)
            ok.append(('mixed kinds raise', False))
        except RulerError:
            ok.append(('mixed kinds raise', True))
    finally:
        shutil.rmtree(d, ignore_errors=True)

    width = max(len(n) for n, _ in ok)
    for name, good in ok:
        print('   %-*s %s' % (width, name, 'ok' if good else 'FAILED'))
    bad = [n for n, g in ok if not g]
    print('\n%d checks, %d failed' % (len(ok), len(bad)))
    return 1 if bad else 0


if __name__ == '__main__':
    import sys
    if len(sys.argv) == 1 or sys.argv[1] == '--selftest':
        sys.exit(_selftest())
    if len(sys.argv) == 2:
        print('%s  %s  [%s]' % (content_md5(sys.argv[1]), sys.argv[1], ruler_for(sys.argv[1])))
    else:
        a, b = sys.argv[1], sys.argv[2]
        print('%s  [%s]' % ('same' if same(a, b) else 'DIFFERS', ruler_for(a)))
