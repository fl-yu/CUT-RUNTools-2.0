#!/usr/bin/python
# Copyright (C) 2019 Qian Zhu
#
# CUT&RUNTools is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.

import sys

f = open(sys.argv[1])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	ct = int(float(ll[-1]))
	if ct==0: continue
	sys.stdout.write("%s\t%s\t%s\t%d\n" % (ll[0], ll[1], ll[2], ct))
f.close()


