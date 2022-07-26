#!/usr/bin/python

import sys

f = open(sys.argv[1])
tt = []
seq = {}
for l in f:
    l = l.rstrip("\n")
    if l.startswith(">"):
        tit = l[1:]
        seq[tit] = ""
        tt.append(tit)
        continue
    seq[tit] += l
f.close()

fw = open(sys.argv[1], "w")
for t in tt:
	c1, c2 = t.split(":")
	c2a,c2b = c2.split("-")
	c2a = int(c2a)+1
	c2b = int(c2b)
	fw.write(">%s:%d-%d\n" % (c1, c2a, c2b))
	fw.write("%s\n" % seq[t])
fw.close()
