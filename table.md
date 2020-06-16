| Hisat2 type | stringtie type| merge type | sensitivity% | precision% | novel loci% | # query transfrags |
|---|---|---|---|---|---|---|
| not strand specific | no parameters | no parameters | 99.9-100 | 62.5-90.5 | 17.5 | 108646 |
| strand specific | --rf -f .1 -c 3 -s 7 -j 2 | no parameters | 100 | 51.2-59.8 | 42.2 | 132608 |
| strand specific | --rf -f .1 -c 3 -s 7 -j 2 | -f 0.1 | 100 | 57.9-64.1 | 42.1 | 112319 |
| strand specific | --rf -f .1 -c 3 -s 7 -j 2 | -f 0.1 -F 1 -T 1 | 100 | 57.9-64.1 | 42.1 | 112319 |
| strand specific | --rf -f .1 -c 3 -s 7 -j 2 | -f 0.1 -F 5 -T 5 | 100 | 74.4-80.2| 23.5| 85417 |
| strand specific | -f .1 -c 3 -s 7 -j 2 | no parameters | 100 | 51.2-59.8 | 42.2 | 132608 |
| strand specific | --rf | no parameters | 100 | 50.2-59.3| 41.4 | 135112 |
| strand specific | -t -c 1.5 -f 0.05 (--conservative no longer a useable flag) | no parameters | 100 | 49.2-58.8 | 40.8 | 137896 |
| strand specific | no parameters | -f 0.1 -F 5 -T 5 | 100 | 74.5-80.3 | 23 | 85265 |
| strand specific | no parameters | no parameters | 100 | 50.2-58.3 | 41.4 | 135112 |
| strand specific | -e | no parameters | 100 | 100 | 0 | 67883 |
