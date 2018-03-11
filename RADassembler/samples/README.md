
An example run:
---
read1: reads containing enzyme cut sites.

read2: randomly sheared reads (paired-end reads).
```
# RADassembler -i read1 -o Assembly_out -s read2 -f gzfastq -P PopMap -t 8 -m 6 -n 6 -M 3 -D 10:100
```
