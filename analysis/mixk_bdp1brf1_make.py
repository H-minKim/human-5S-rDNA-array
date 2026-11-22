# -*- coding: utf-8 -*-
# Writes mixk_bdp1brf1.py : mixk_run.py with its tracks replaced by K562 BDP1 and BRF1 ChIP-seq, the same marker lists
#   and the same gate, writing K562_QUALGATE_MIX_BDP1BRF1.tsv with scratch files of its own. Run on the compute node.
import io, os
W = os.path.expanduser("~/p11_encode")
s = io.open(os.path.join(W, "mixk_run.py"), encoding="utf-8").read()
a = s.index("JOBS = ["); b = s.index("]\n", a) + 2
s = s[:a] + ('JOBS = [\n    ("BDP1",    [RG % x for x in ("ENCFF000YGF", "ENCFF000YGH", "ENCFF000YGI")]),\n'
             '    ("BRF1",    [RG % x for x in ("ENCFF000YGP", "ENCFF000YGR", "ENCFF000YGS")]),\n]\n') + s[b:]
assert s.count('"K562_QUALGATE_MIX.tsv"') == 1
s = s.replace('"K562_QUALGATE_MIX.tsv"', '"K562_QUALGATE_MIX_BDP1BRF1.tsv"').replace("jbmix_", "jbmixBB_")
io.open(os.path.join(W, "mixk_bdp1brf1.py"), "w", encoding="utf-8", newline="\n").write(s)
compile(s, "x", "exec"); print("mixk_bdp1brf1.py written")
