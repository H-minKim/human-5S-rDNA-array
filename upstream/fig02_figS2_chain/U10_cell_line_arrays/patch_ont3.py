# -*- coding: utf-8 -*-
# Transcribed from the run record: the code ran as a shell heredoc or a one-line command, so no script file existed on disk.
# Patch applied on the node to k562_ont3.sh ( line 60330 written, applied at line 60341, 2026-09-23 07:11 UTC).
import os
p = os.path.expanduser("~/p11_encode/k562_ont3.sh")
s = open(p).read()
old = '  if [ "$got" != "$md5" ]; then'
new = '  if [ "$md5" = "NOMD5" ]; then log "$acc md5 모름 - 크기만 확인"; elif [ "$got" != "$md5" ]; then'
if new.split(";")[0] in s:
    print("이미 패치됨")
else:
    assert old in s, "패턴 못 찾음"
    open(p, "w").write(s.replace(old, new, 1))
    print("패치 완료")
