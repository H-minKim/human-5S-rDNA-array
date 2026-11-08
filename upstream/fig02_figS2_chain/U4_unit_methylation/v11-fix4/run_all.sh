#!/bin/bash
# v11-fix4: reproduce only HG008-T from the local BAM (CPGONLY=0), and CpG only (CPGONLY=1)
N="${RDNA_V11_FIX4_DIR:?set RDNA_V11_FIX4_DIR}"; C="$N/01_code/local"; cd "$C"
for k in ao jx; do
  CPGONLY=0 PYTHONIOENCODING=utf-8 python meth_units_${k}_v11_CPGONLY.py "$N/09_tests/repro_CPGONLY0/out_$k" LIST_HG008T.txt > "$N/09_tests/repro_CPGONLY0/$k.log" 2>&1 &
  CPGONLY=1 PYTHONIOENCODING=utf-8 python meth_units_${k}_v11_CPGONLY.py "$N/09_tests/CPGONLY1/out_$k" LIST_HG008T.txt > "$N/09_tests/CPGONLY1/$k.log" 2>&1 &
done
wait; echo DONE
