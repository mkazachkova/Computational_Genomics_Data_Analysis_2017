"""This one takes ~1 hour to run!"""

import sys
import urllib2
from bs4 import BeautifulSoup
import numpy as np


snp = open('data_input/GSE36138_series_matrix.txt')
for i in range(69):
	snp.readline() #trash

headers_snp = snp.readline().strip().split('\t')
headers_snp.pop(0)

#take in data for snp array 
count = 0
data_snp = np.empty((0,947), float)
for line in snp:
	temp = line.strip().split('\t')
	if len(temp) < 100:
		continue
	temp.pop(0)
	temp = [temp]
	data_snp = np.append(data_snp, temp, axis=0)
	print count
	count+=1


data_snp = data_snp.transpose()
np.save('temp_output/snp_array.npy', data_snp) #will be saved in format cell lines x samples



gene = open('data_input/GSE36133_series_matrix.txt')
for i in range(69):
	gene.readline() #trash 

headers_gene = gene.readline().strip().split('\t')
headers_gene.pop(0)
print headers_gene

#take in data for gene expression
data_gene = np.empty((0,917), float)


count = 0
for line in gene:
	temp = line.strip().split('\t')
	if len(temp) < 100:
		continue
	temp.pop(0)
	temp = [temp]
	data_gene = np.append(data_gene, temp, axis=0)
	print count
	count+=1


data_gene = data_gene.transpose()
np.save('temp_output/gene_array.npy', data_gene) #will be saved in format cell lines x samples
