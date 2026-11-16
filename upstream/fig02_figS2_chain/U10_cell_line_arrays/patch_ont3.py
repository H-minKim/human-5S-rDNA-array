# -*- coding: utf-8 -*-
# Patch applied on the node to k562_ont3.sh (tolerates a missing md5).
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
