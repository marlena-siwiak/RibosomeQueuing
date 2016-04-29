#!/usr/bin/env python
# -*- coding: utf-8 -*-

VERSION = 3.0

from fasta import *
from os import mkdir, system
import sys

if len(sys.argv) != 6:
	print "\nUsage:\t%s\tfile.fasta\ttRNA_time.tsv\topt_codons\tmutated_seq\tmutations_rep\n" % sys.argv[0]
	sys.exit(1)


ORFS = sys.argv[1]
TRNA_CZAS = sys.argv[2]
OPT_CODONS = sys.argv[3]
FOLDER_SEQMUT = sys.argv[4]
FOLDER_RAPMUT = sys.argv[5]


##### Helping objects.

class pozycjeryb(list):
	"""Shows ribosome positions in the sequence."""
	def __init__(self):
		list.__init__(self)
	def appendp(self, poz):
		for p in range(poz-5, poz+6):
			if p >= 0:
				self.append(p)
	def removep(self, poz):
		for p in range(poz-5, poz+6):
			if p >= 0:
				self.remove(p)
	def countp(self, poz):
		"""Returns 1 if at least one of the flanking positions is covered by a ribosoe."""
		for p in range(poz-5, poz+6):
			if p >= 0:
				if self.count(p) != 0:
					return 1
		return 0
	
class polisomy(dict):
	"""Contains information on ribosomes' positions."""
	def __init__(self, seq):
		dict.__init__(self)
		self.pozycje = pozycjeryb()				
		self.seq = seq						
		self.czasy = {}						
		self.evilcod = {}					
	def __get_wait(self, poz):
		"""Returns the required time of waiting in position poz. If poz is outside the coding sequence or on a stop codon, returns 0."""
		try:
			kodon = self.seq[poz]
		except IndexError:
			wait = 0			
		else:
			if trnadict.has_key(kodon):
				wait = int(round(trnadict[kodon], 0))
			else:
				wait = 0
		return wait
	def __move(self, rybosom):
		"""If possible, moves ribosome one step forward, poz refers to the codon currently being translated."""
		poz = self[rybosom][0]				
		if self.pozycje.count(poz+6) == 0:			
			wait = self.__get_wait(poz+1)
			if wait == 0:					
			# if the sequence is over, remove ribosome and save result.
				wynik = self.pop(rybosom)
				self.czasy[rybosom] = wynik[3]	
				self.pozycje.removep(poz)
			else:
			# otherwise, move the ribosome.	
				self[rybosom][0] += 1
				self[rybosom][1] = 0
				self[rybosom][2] = wait
				self.pozycje.removep(poz)
				self.pozycje.appendp(poz+1)
		elif self.evilcod.has_key(poz + 11) == False:
			self.evilcod[poz + 11] = 1
		elif self.evilcod.has_key(poz + 11):
			self.evilcod[poz + 11] += 1
	def add_rib(self, numer):
		"""Adds new ribosome to the mRNA sequence."""
		name = "rybosom" + str(numer)
		wait = self.__get_wait(0)
		if self.pozycje.countp(0) == 0:				
			self[name] = [0, 0, wait, 0]
			self.pozycje.appendp(0)
		else:
			raise IndexError
	def move_rib(self):
		"""If possible moves ribosomes one position forward."""
		keys = self.keys()
		keys.sort(reverse=True)					
		for r in keys:
			self[r][1] = self[r][1]+1			
			self[r][3] = self[r][3]+1			
			if self[r][1] >= self[r][2]:			
				self.__move(r)

##### Helping functions.

def make_trnadict(plik):
	"""Prepares a dict with translation times of codons from a given file."""
	trnafile = open(plik, "rU")
	trnadict = {}
	for line in trnafile:
		if not line.startswith("kodon"):
			line = line.strip()
			line = line.split("\t")
			trnadict[line[0]] = float(line[1].strip())
	trnafile.close()
	return trnadict

def elongation_time(seq):
	"""Calculates elongation time of a sequence."""
	czas = 0
	for s in seq:
		if trnadict.has_key(s):
			czas = czas + round(trnadict[s], 0)
	czas = int(round(czas, 0))
	return czas


def translacja(seq, tini):
	"""Simulates the translation process. Returns the minimal possible time between subsequent initiations of ribosomes for which no ribosome queueing occurs. Returns the difference betwwen the elongation time of the first and second ribosome, if the ribosomes initiate as with the highest possible frequency."""
	polisom = polisomy(seq)
	nr = 1								# ribosome number
	ni = 0								# time to the next initiation
	while len(polisom.czasy) != 2:					# its enough to simulate for two ribosomes only
		if ni == 0 and nr <= 2:					
			try:			
				polisom.add_rib(nr)
			except IndexError:				# the 5' end is blocked
				print "To sie nigdy nie powinno wydrukowac. Jesli sie wydrukowalo, trzeba zweryfikowac rozumowanie."
				sys.exit(1)
			else:
				nr += 1
				ni = tini + 1				
		polisom.move_rib()
		ni = ni - 1
	return polisom.evilcod

def count_ie(evilcod):
	"""Returns the minimal possible Ie or 0 if there is no queueing."""
	suma = 0
	for v in evilcod.values():
		suma = suma + v
	return suma

def codque(evilcod):
	"""Returns the number of the first codon that causes queuing, or -1 if there is no queuing."""
	if evilcod != {}:
		cod = evilcod.keys()
		cod.sort()
		return cod[0]
	else:
		return -1

def codmut(seq, nr_cod):
	"""Returns the number of codon to mutate. Returns -1 if there is no such codon."""
	checklist = list(seq[0:nr_cod+1][::-1])
	for c in checklist:
		if optdict.has_key(c):
			return nr_cod
		else:
			nr_cod = nr_cod - 1
	return -1

def mutate(seq, nr_cod):
	"""Mutates the codon of a given nr_cod. Returns the mutated sequence."""
	cod = seq[nr_cod]
	newseq = list(seq[0:nr_cod])
	newseq.append(optdict[cod])
	newseq += seq[nr_cod+1:]
	return newseq
	


##### Main.

# Create the fasta object for sequences.
mrnadict = makefasta(ORFS)
mrnadict.reducekey()
mrnadict.touracyl()
mrnadict.tomrna()

# Create a dict with codons' translation times.
trnadict = make_trnadict(TRNA_CZAS)

# Create a dict for mutating codons and a dict of translation code.
optdict = {}
aadict = {}
f = open(OPT_CODONS, "rU")
for line in f:
	if not line.startswith("aa"):
		line = line.strip()
		line = line.split("\t")
		kodony = line[2].split(",")
		for k in kodony:
			if not k == line[1]:
				optdict[k] = line[1]
			aadict[k] = line[0]
f.close()
# Create output folders for fasta files and report files for each sequence.
for folder in (FOLDER_SEQMUT, FOLDER_RAPMUT):
	try:
		mkdir(folder)				
	except OSError:
		system("rm -r " + folder)
		mkdir(folder)


print "gene\tnr_mut\tproc_mut\tIe_org\tIq_org\tImin_org\tIe_ost\tIq_ost\tImin_ost"

for gen in mrnadict.keys():
	seq    = mrnadict[gen]
	mutseq = fasta()			
	rapfile  = open(FOLDER_RAPMUT+"/"+gen+"_rapmut.tsv", "w")		
	rap_head = ["seq_name", "nr_codque", "typ_codque", "aa_codque", "nr_codmut", "typ_codmut", "aa_codmut", "Ie", "Iq", "Imin"]
	rapfile.write("\t".join(rap_head) + "\n")
	nr_mut = 0				
	seq_name   = gen + "_mut_" + str(nr_mut)
	Iq      = elongation_time(seq[0:11])
	evilcod = translacja(seq, Iq)
	Ie      = count_ie(evilcod)
	Imin    = Iq + Ie

	nr_codque  = codque(evilcod)
	nr_codmut  = codmut(seq, nr_codque)

	Ie_org   = Ie
	Iq_org   = Iq
	Imin_org = Imin
	
	while nr_codque != -1 and nr_codmut != -1:
		typ_codque = seq[nr_codque]
		typ_codmut = seq[nr_codmut]
		aa_codque  = aadict[typ_codque]
		aa_codmut  = aadict[typ_codmut]
		seq_name   = gen + "_mut_" + str(nr_mut)
		mutseq[seq_name] = seq
		rap_line = [seq_name, str(nr_codque + 1), typ_codque, aa_codque, str(nr_codmut + 1), typ_codmut, aa_codmut, str(Ie), str(Iq), str(Imin)]
		rapfile.write("\t".join(rap_line) + "\n")
		seq = mutate(seq, nr_codmut)
		nr_mut += 1
		evilcod = translacja(seq, Iq)
		Iq = elongation_time(seq[0:11])
		Ie = count_ie(evilcod)
		Imin = Iq + Ie
		nr_codque  = codque(evilcod)
		nr_codmut  = codmut(seq, nr_codque)

	seq_name   = gen + "_ost"
		
	if nr_codque == -1:
		rap_line = [seq_name, "NA", "NA", "NA", "NA", "NA", "NA", str(Ie), str(Iq), str(Imin)]
	elif nr_codmut == -1:
		typ_codque = seq[nr_codque]
		aa_codque  = aadict[typ_codque]
		rap_line = [seq_name, str(nr_codque + 1), typ_codque, aa_codque, "NA", "NA", "NA", str(Ie), str(Iq), str(Imin)]
		
	mutseq[seq_name] = seq
	mutseq.printfasta(FOLDER_SEQMUT+"/"+gen+"_mut.fa")
	
	rapfile.write("\t".join(rap_line) + "\n")
	rapfile.close()

	proc_mut = round(float(nr_mut)/len(seq)*100, 0)
	out_line = [gen, str(nr_mut), str(proc_mut), str(Ie_org), str(Iq_org), str(Imin_org), str(Ie), str(Iq), str(Imin)]
	print "\t".join(out_line)

