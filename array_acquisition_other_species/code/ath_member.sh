#!/bin/bash
# tar --to-command 가 멤버마다 부른다 (표준입력 = 멤버 내용 · TAR_FILENAME = 멤버 이름)
#   pod5 멤버 : 임시로 받아 IDS 에 있는 read 만 sub/ 에 남기고 원본은 지운다 -> 디스크에는 고른 read 만
#   members.tsv : 이름 · 크기 · pod5 filter 종료코드 · 멤버 안 read 수 · 남긴 read 수
set -u
H=~/xsp/ath
mkdir -p /tmp/athm "$H/sub"
b=$(basename "$TAR_FILENAME")
case "$b" in
  *.pod5)
    t=/tmp/athm/$b
    cat > "$t"
    sz=$(stat -c %s "$t")
    ~/xsp/venv/bin/pod5 filter "$t" --ids "$IDS" --output "$H/sub/$b" --missing-ok --force-overwrite > "$H/sub/$b.log" 2>&1
    rc=$?
    nin=$(~/xsp/venv/bin/python -c "import pod5,sys; print(len(pod5.Reader(sys.argv[1]).read_ids))" "$t" 2> /dev/null)
    nout=0
    [ -s "$H/sub/$b" ] && nout=$(~/xsp/venv/bin/python -c "import pod5,sys; print(len(pod5.Reader(sys.argv[1]).read_ids))" "$H/sub/$b" 2> /dev/null)
    printf '%s\t%s\t%s\t%s\t%s\n' "$b" "$sz" "$rc" "${nin:-NA}" "${nout:-NA}" >> "$H/members.tsv"
    rm -f "$t"
    ;;
  *)
    cat > /dev/null
    printf '%s\tskip\n' "$b" >> "$H/members.tsv"
    ;;
esac
