#!/usr/bin/python
# Copyright (C) 2019 Qian Zhu
#
# CUT&RUNTools is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.

import sys
import re

ch_size = {}
f = open(sys.argv[1])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	ch_size[ll[0]] = int(ll[1])
f.close()

f = open(sys.argv[2])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	st = int(ll[1])
	ed = int(ll[2])
	chrom = ll[0]
	if chrom==".": continue
	if st>ch_size[chrom] or ed>ch_size[chrom]:
		continue

	chrom2 = ll[3]
	st2 = int(ll[4])
	ed2 = int(ll[5])
	if chrom2==".": continue
	if st2>ch_size[chrom] or ed2>ch_size[chrom]:
		continue

	if chrom!=chrom2:
		continue

	sys.stdout.write(l + "\n")
f.close()

