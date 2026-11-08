# -*- coding: utf-8 -*-
# 규격 러너 - run_panel_split2.py 와 같되 paper_spec_hook 을 먼저 끼운다.
# 정본 코드는 안 건드린다 : 사본을 모래상자에 만들고 출력 경로만 바꾼다.
# 쓰기 : python run_spec_split.py <K> <코드.py> [코드.py ...]
import io
import os
import re
import subprocess
import sys

NL = chr(10)
SP = os.path.dirname(os.path.abspath(__file__))
# ㉤ : 이름만 주면 저장소 안에서 찾는다.
# It used to point at the code folder (01_code/data_generation_analysis)
# and silently read the original rather than the repository copy ("the file name is the interface").
REPO = os.path.dirname(os.path.dirname(SP))
CODE = [os.path.join(REPO, x) for x in
        ("fig04", "figS4", "figS5", "figS6",
         os.path.join("common", "fig04_figS4"), "analysis")]
SB = os.path.join(SP, os.environ.get("PANEL_OUT", "panelsplit_spec"))
BB = os.path.join(SB, "_panels")
for d in (SB, BB):
    os.makedirs(d, exist_ok=True)

K = sys.argv[1]
ONLY = sys.argv[2:]
if not ONLY:
    raise SystemExit("쓰기 : run_spec_split.py <K> <코드.py> ...")
env = dict(os.environ)
env["SPEC_K"] = K
env["SPEC_FLOOR_LOG"] = os.path.join(SB, "_floor.tsv")
env["PANEL_SPLIT_DIR"] = BB
env["PYTHONIOENCODING"] = "utf-8"

# 훅은 모래상자 안에 있어야 사본이 찾는다 (sys.path 에 SP 를 먼저 넣는다)
for fn in ONLY:
    src = fn if os.path.isabs(fn) else None
    if src is None:
        for _d in CODE:
            if os.path.exists(os.path.join(_d, fn)):
                src = os.path.join(_d, fn)
                break
    fn = os.path.basename(fn)
    if src is None or not os.path.exists(src):
        print("%-38s 저장소에 없음 (절대경로로 주십시오)" % fn)
        continue
    s = io.open(src, encoding="utf8").read()
    out = os.path.join(SB, fn[:-3])
    os.makedirs(out, exist_ok=True)

    def repl(m, _o=out):
        return "%s = r\"%s\"" % (m.group(1), _o)

    s2 = re.sub(r"^(OUTG|OUTD|OUTT|OUT|CLIP)\s*=\s*.+$", repl, s, flags=re.M)
    head = ("import sys as _s" + NL
            + "_s.path.insert(0, r\"%s\")" % SP + NL
            + "import os as _o" + NL
            # ㉥ : 사본은 모래상자에서 돌아 __file__ 이 원본 폴더를
            # 가리키지 않는다 ⇒ 원본 폴더를 물려준다 (입력을 그 옆에서 찾는 글이 있다).
            + "_o.environ[\"SPEC_SRC_DIR\"] = r\"%s\"" % os.path.dirname(src) + NL
            + "_o.environ[\"PANEL_SPLIT_DIR\"] = r\"%s\"" % BB + NL
            + "import matplotlib as _m; _m.use(\"Agg\")" + NL
            + "import paper_spec_hook  # noqa  먼저" + NL
            + "import panel_split_hook2  # noqa  그 다음" + NL
            + "_m.rcParams[\"svg.fonttype\"] = \"none\"; _m.rcParams[\"pdf.fonttype\"] = 42" + NL)
    p = os.path.join(SB, fn)
    io.open(p, "w", encoding="utf8", newline=NL).write(head + s2)
    r = subprocess.run([sys.executable, "-X", "utf8", p], cwd=SB, env=env,
                       capture_output=True, text=True, encoding="utf8",
                       errors="replace", timeout=5400)
    nb = len([x for x in os.listdir(BB) if x.endswith(".panels.json")])
    print("%-38s %-6s 누적 %2d  %s"
          % (fn, "OK" if r.returncode == 0 else "실패", nb,
             (r.stderr.strip().split(NL)[-1][:90] if r.returncode else "")))
    sys.stdout.flush()
print("패널 : %s" % BB)
