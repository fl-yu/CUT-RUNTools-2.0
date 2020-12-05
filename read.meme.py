#!/usr/bin/python
import sys
import os
import re

def read_summary(n):
	f = open(n)
	h = f.readline().rstrip("\n").split("\t")
	alist = []
	for l in f:
		l = l.rstrip("\n")
		if l.startswith("#"):
			continue
		if l=="":
			continue
		ll = l.split("\t")
		cx = dict(zip(h, ll))
		alist.append(cx)
	f.close()
	return alist

def read_meme(n):
	f = open(n)
	st = 1
	matrix = {}
	while True:
		l = f.readline()
		if l=="": break
		l = l.rstrip("\n")
		if l.find("MEME-%d position-specific probability matrix" % st)!=-1:
			matrix.setdefault("MEME-%d" % st, [])
			f.readline()
			h = f.readline().rstrip("\n")
			tmp = h.find(" w= ")
			tmpstr = h[tmp+4:].split()
			length = int(tmpstr[0])
			matrix["MEME-%d" % st].append(h)
			for i in range(length):
				l1 = f.readline().rstrip("\n")
				matrix["MEME-%d" % st].append(l1)
			#print st, length
			st+=1
	f.close()
	return matrix

def read_dreme(n):
	f = open(n)
	st = 1
	matrix = {}
	while True:
		l = f.readline()
		if l=="": break
		l = l.rstrip("\n")
		if l.find("DREME-%d" % st)!=-1:
			matrix.setdefault("DREME-%d" % st, [])
			f.readline()
			while True:
				l1 = f.readline().rstrip("\n")
				if l1.startswith("# "):
					continue
				if l1=="": break
			h = f.readline().rstrip("\n")
			tmp = h.find(" w= ")
			tmpstr = h[tmp+4:].split()
			length = int(tmpstr[0])
			matrix["DREME-%d" % st].append(h)
			for i in range(length):
				l1 = f.readline().rstrip("\n")
				matrix["DREME-%d" % st].append(l1)
			st+=1
	f.close()
	return matrix
	
def print_motif_meme(motif, matrix):
	s = []
	s.append("MEME version 4")
	s.append("")
	s.append("ALPHABET= ACGT")
	s.append("")
	s.append("strands: + -")
	s.append("")
	s.append("Background letter frequencies (from uniform background):")
	s.append("A 0.25000 C 0.25000 G 0.25000 T 0.25000")
	s.append("")
	s.append("MOTIF 1 " + motif)
	s.append("")
	s.extend(matrix)
	return "\n".join(s)
	
if __name__=="__main__":
	this_dir = sys.argv[1]
	ss = read_summary(this_dir + "/summary.tsv")
	meme_matrices = read_meme(this_dir + "/meme_out/meme.txt")
	dreme_matrices = read_dreme(this_dir + "/dreme_out/dreme.txt")
	outdir = this_dir + "/motifs"
	
	if not os.path.isdir(outdir):
		sys.stderr.write("Creating directory in %s...\n" % outdir)
		os.mkdir(outdir)

	sys.stderr.write("Motifs created in directory %s\n" % outdir)
	for ind,cx in enumerate(ss):
		this_id = cx["ALT_ID"]
		this_src = cx["MOTIF_SOURCE"]
		this_eval = cx["E-VALUE"]
		this_motif = cx["CONSENSUS"]
		filename = outdir + "/" + this_id + "." + this_motif + ".meme"
		if this_src=="MEME":
			fw = open(filename, "w")
			this_str = print_motif_meme(this_motif, meme_matrices[this_id])
			fw.write(this_str)
			fw.close()
		elif this_src=="DREME":
			fw = open(filename, "w")
			this_str = print_motif_meme(this_motif, dreme_matrices[this_id])
			fw.write(this_str)
			fw.close()
		#print this_id, this_src, this_eval, this_motif
		#if this_src=="MEME":
