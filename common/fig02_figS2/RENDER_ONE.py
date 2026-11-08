# -*- coding: utf-8 -*-
r"""Redraw one figure in Korean and English (in the mirror only). Usage: python RENDER_ONE.py FIG_XXX.py tag"""
import io, os, re, sys, time, shutil, subprocess
sys.stdout.reconfigure(encoding="utf-8")
HERE = os.path.dirname(os.path.abspath(__file__))
P = (__import__("os").environ["RDNA_DATA_ROOT"] + r"/for_paper/00_paper_figures_distal_layout")
CODE = os.path.join(P, "code")
REPRO = os.path.join(HERE, "REPRO")
FIX = os.path.join(HERE, "_fix")
CODEF = sys.argv[1]
TAG = sys.argv[2]
FT = os.path.join(REPRO, "figtab_h1f")
DEFAULT = 110
BUD = {}
for ln in io.open(os.path.join(REPRO, "EN_WRAP_MAP.tsv"), encoding="utf8"):
    if not ln.startswith("#") and "\t" in ln:
        f = ln.rstrip("\n").split("\t"); BUD[f[0]] = int(f[1])


def prep(tmp):
    assert not tmp.startswith(("F:", "f:"))
    if os.path.exists(tmp): shutil.rmtree(tmp)
    os.makedirs(tmp)
    for f in os.listdir(CODE):
        if f.endswith(".py"): shutil.copyfile(os.path.join(CODE, f), os.path.join(tmp, f))
    for f in os.listdir(FIX):
        if f.endswith(".py"): shutil.copyfile(os.path.join(FIX, f), os.path.join(tmp, f))
    os.makedirs(os.path.join(tmp, "REPRO"), exist_ok=True)
    for f in os.listdir(FT):
        if f.endswith(".tsv"): shutil.copyfile(os.path.join(FT, f), os.path.join(tmp, "REPRO", f))
    for sub in ("_vt2", "_vt3"):
        s = os.path.join(REPRO, sub)
        if os.path.isdir(s):
            d = os.path.join(tmp, "REPRO", sub); os.makedirs(d, exist_ok=True)
            for f in os.listdir(s):
                if f.endswith(".tsv"): shutil.copyfile(os.path.join(s, f), os.path.join(d, f))


MISS = os.path.join(REPRO, "EN_MISS_%s.tsv" % TAG)
if os.path.exists(MISS): os.remove(MISS)
_HK = os.environ.get("RENDER_HOOKS", "hook_sym,hook_en_sym").split(",")   # the MAJREF family uses hook_tieronly,hook_tieronly_en
for lang, hook, budget in (("ko", _HK[0], False), ("en", _HK[1], True)):
    tmp = os.path.join(HERE, "_t1_%s_%s" % (TAG, lang))
    prep(tmp)
    out = os.path.join(REPRO, "_r%s%s" % (TAG, lang))
    assert not out.startswith(("F:", "f:"))
    if os.path.exists(out): shutil.rmtree(out)
    os.makedirs(out)
    e = {**os.environ, "PYTHONIOENCODING": "utf-8", "PYTHONHASHSEED": "0",
         "FIGTAB": FT, "FIGOUT": out, "OUTDIR": out,
         "V11_ROOT": os.path.join(REPRO, "mirror_h1f"),
         "REPRO_MIRROR": os.path.join(REPRO, "mirror_h1f"),
         "V11_TIER": "h1", "PTABDIR": os.path.join(REPRO, "ptab_h1"),
         "TRIOTAB": os.path.join(REPRO, "out_trio_maj"),
         "PYTHONPATH": os.path.join(HERE, hook) + os.pathsep + os.environ.get("PYTHONPATH", "")}
    e.pop("EN_WIDTH_LOG", None)
    if budget:
        e["EN_WRAP_MAXLINE"] = str(BUD.get(CODEF, DEFAULT))
        e["EN_DICT"] = os.path.join(HERE, "EN_DICT.tsv")
        e["EN_MISS"] = MISS
    t0 = time.time()
    cmd = [sys.executable, os.path.join(tmp, CODEF)]
    if os.environ.get("ART_DUMP"):          # write what was drawn as JSON (for content comparison), the figure output is the same
        cmd = [sys.executable, os.path.join(tmp, "_dumpwrap.py"), os.path.join(tmp, CODEF), os.path.join(out, "_art.json")]
    if os.environ.get("RETARGET_SPEC"):     # supplementary cell figure, run the original code as is and move only one panel into the cell frame
        cmd = [sys.executable, os.path.join(tmp, "_retarget.py"), os.path.join(tmp, CODEF),
               os.environ["RETARGET_SPEC"].replace("{L}", lang), out]
    if os.environ.get("TRACE_OPEN"):        # record the files actually read (for the code repository list, the figure output is the same)
        os.makedirs(os.environ["TRACE_OPEN"], exist_ok=True)
        cmd = [sys.executable, os.path.join(tmp, "_tracewrap.py"), os.path.join(os.environ["TRACE_OPEN"], "%s__%s__%s.tsv" % (TAG, CODEF, lang))] + cmd[1:]
    if os.environ.get("NOHOOK_FIGLANG"):    # code that picks the language with FIG_LANG without a hook, like FIG_COSHARE
        e["PYTHONPATH"] = os.environ.get("PYTHONPATH", ""); e["FIG_LANG"] = lang; e.pop("EN_DICT", None)
    r = subprocess.run(cmd, cwd=tmp, env=e,
                       capture_output=True, text=True, encoding="utf8", errors="replace", timeout=3600)
    ok = r.returncode == 0
    print("   %-3s %-22s %s %5.1fs" % (lang, CODEF, "ok" if ok else "실패", time.time() - t0))
    _all = (r.stdout or "").strip().split("\n")
    for ln in [x for x in _all[:-6] if "닿음 검사" in x or "칸 밖 글자" in x or "가장 가까운 짝" in x or "본 글자" in x] + _all[-6:]:   # all check lines
        if ln.strip(): print("      | %s" % ln[:260])
    if not ok:
        print("      %s" % (r.stderr or "").strip().split("\n")[-1][:220])
        sys.exit(1)
    print("      산출 : %s" % ", ".join(sorted(os.listdir(out))[:6]))
HAN = re.compile("[가-힣]")
for f in sorted(os.listdir(os.path.join(REPRO, "_r%sen" % TAG))):
    if f.endswith(".svg"):
        s = io.open(os.path.join(REPRO, "_r%sen" % TAG, f), encoding="utf8").read()
        han = {x.strip() for x in re.findall(r">([^<]*)<", s, re.S) if HAN.search(x)}
        print("   영어판 %-34s 한글 %d %s" % (f, len(han), sorted(han)[:4]))
if os.path.exists(MISS):
    uniq = {l.rstrip("\n").split("\t")[-1] for l in io.open(MISS, encoding="utf8")}
    print("   번역 안 된 한글 %d : %s" % (len(uniq), sorted(uniq)[:8]))
