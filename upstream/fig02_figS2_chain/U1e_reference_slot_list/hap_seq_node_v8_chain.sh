#!/bin/bash
# hap_seq_node_v8.tsv (616 rows: sample, slot, reference file) is the reference list read by the unit map (U2).
# Step 0, seqassign.py (hap_seq.tsv), is analysis/refset_step3_assign_slots.py
# Steps 1-3 below are the lines that write hap_seq_node*.tsv. They need the node files.
# Machine paths are given by RDNA_DATA_ROOT / RDNA_WORK_DIR (no default).
# step 1 (local): seqassign output hap_seq.tsv -> node list
tail -n +2 hap_seq.tsv | awk -F'\t' '$7!="-"{print $1"\t"$2"\t"$7}' > hap_seq_node.tsv && wc -l < hap_seq_node.tsv && head -3 hap_seq_node.tsv

# step 2 (local, G = folder holding hap_seq_node.tsv): six slots to the read-consensus reference -> v7
G="$RDNA_WORK_DIR/g8v2"
head -2 "$G/hap_seq_node.tsv"; wc -l < "$G/hap_seq_node.tsv"
awk -F'\t' 'BEGIN{OFS="\t"} {k=$1"_"$2} k=="HG02683_1"||k=="NA19090_1"||k=="HG03126_1"||k=="HG03366_1"||k=="NA19120_1"||k=="HG03198_1"{$3="cons2seq__58__c2seq_"k".fa"} {print}' "$G/hap_seq_node.tsv" > "$G/hap_seq_node_v7.tsv"
diff "$G/hap_seq_node.tsv" "$G/hap_seq_node_v7.tsv" | grep '^>'

# step 3 (on the node): hg002 slot 2 -> HG002 Q100 paternal array -> v8
mkdir -p ~/cmp4d_v8
for f in ~/cmp4d_v7/*.fa; do ln -sf "$f" ~/cmp4d_v8/$(basename "$f"); done
cp ~/hg002asm/pat_array.fa ~/cmp4d_v8/hg002asm__58__pat_array.fa
echo "=== cmp4d_v8 파일수"; ls ~/cmp4d_v8 | wc -l
python3 - <<PY
import io
NL=chr(10); TAB=chr(9)
src="$NODE_HOME/hap_seq_node_v7.tsv"; dst="$NODE_HOME/hap_seq_node_v8.tsv"
out=[]; n=0
for L in io.open(src,encoding="utf8"):
    f=L.rstrip(NL).split(TAB)
    if len(f)>=3 and f[0]=="hg002" and f[1]=="2":
        f[2]="hg002asm__58__pat_array.fa"; n+=1
    out.append(TAB.join(f))
io.open(dst,"w",encoding="utf8",newline=NL).write(NL.join(out)+NL)
print("hg002 슬롯2 교체", n, "줄 · 전체", len(out), "줄")
PY
grep -P "^hg002\t" ~/hap_seq_node_v8.tsv
