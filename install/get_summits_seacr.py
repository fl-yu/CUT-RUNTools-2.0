#!/usr/bin/python
# Copyright (C) 2019 Qian Zhu
#
# CUT&RUNTools is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; version 2 of the License.
#
# CUT&RUNTools is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# A copy of the GNU General Public License has been distributed along with CUT&RUNTools and is found in LICENSE.md.

import sys
import os
import re

f = open(sys.argv[1])
for l in f:
	l = l.rstrip("\n")
	ll = l.split("\t")
	coord = ll[-1]
	c1 = coord.split(":")[0]
	cx = coord.split(":")[1]
	t1 = int(cx.split("-")[0])
	if cx.split("-")[1]=="":
		t2 = t1
	else:
		t2 = int(cx.split("-")[1])
	
	mid = (t1+t2)/2
	sys.stdout.write("%s\t%d\t%d\n" % (c1, mid, mid+1))
f.close()
