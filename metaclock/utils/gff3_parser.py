#!/usr/bin/env python
from collections import namedtuple


def attribute_dict(attri):
	dict_={}
	for a in attri.split(';'):
		k=a.split("=")[0]
		v=a.split("=")[1]
		dict_[k]=v
	return dict_	


class Parser(object):

	def __init__(self, gff3):
		self.gff3 =(i.rstrip() for i in open(gff3).read().split("##")[2].rstrip().split('\n')[1:])

	def element(self):
		Element=namedtuple('Element','seqid source type start end score strand phase attributes')
		for I in self.gff3:
			i = I.split('\t')	
			yield Element(seqid = i[0], source = i[1], type = i[2], start = i[3], end = i[4],
				         score = i[5], strand = i[6], phase = i[7], attributes = attribute_dict(i[8]))
		


