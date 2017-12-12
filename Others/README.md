A modified version of tsv2bam from [Stacks2](http://catchenlab.life.illinois.edu/stacks/stacks_v2.php)
=
modifications: support more fastaq header formats

* type I: end with 1 or 2 seperated by a delimator (|, / or _), e.g. <b>@4_1101_13393_18/1</b>

* type II: seperated by space, containing index informations, e.g. <b>@HWI-ST1276:71:C1162ACXX:1:1101:1208:2458 1:N:0:CGATGT</b>

How to install
---
```
tar -xvf stacks-2.0Beta6-modified.tar.gz
cd stacks-2.0Beta6 && make -f Tsv2bam.make
```

Copy the new compiled tsv2bam, then run!
