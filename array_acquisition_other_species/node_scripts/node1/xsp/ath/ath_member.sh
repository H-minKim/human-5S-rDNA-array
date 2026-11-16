#!/bin/bash
# tar --to-command calls this once per member: stdin is the member's content and TAR_FILENAME is its name.
#   For a pod5 member: take it temporarily, keep only the reads listed in IDS under sub/, then delete the original, so only the selected reads reach disk.
#   members.tsv: name, size, pod5 filter exit code, reads in the member, reads kept.
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
