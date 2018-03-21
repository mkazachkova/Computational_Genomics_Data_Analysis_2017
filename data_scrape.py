import numpy as np
import urllib2
from bs4 import BeautifulSoup

result_snp = np.load('data_output/snp_array.npy')
#print result_snp.shape
result_gene = np.load('data_output/gene_array.npy')
#print result_gene.shape



snp = open('data_input/GSE36138_series_matrix.txt')
for i in range(69):
	snp.readline() #trash

headers_snp = snp.readline().strip().split('\t')
headers_snp.pop(0)
print len(headers_snp)

gene = open('data_input/GSE36133_series_matrix.txt')
for i in range(69):
	gene.readline() #trash 

headers_gene = gene.readline().strip().split('\t')
headers_gene.pop(0)
print len(headers_gene)


cell_lines = open('data_input/CCLE_NP24.2009_Drug_data_2015.02.24.csv')
cell_lines.readline() #trash
lines = []
ics = []
for line in cell_lines:
	line = line.strip().split(',')
	lines.append(line[0])
	ics.append(line[len(line) - 3])


final_snp_data = np.empty((0,22419), float)
final_gene_data = np.empty((0,18926), float)
final_cell_lines = []
final_ics = []



j = 0
for i in range(len(headers_snp)):
	temp = headers_snp[i][1:(len(headers_snp[i])-1)]
	#break
	link1 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + temp
	page = urllib2.urlopen(link1)
	soup = BeautifulSoup(page)
	t =  soup.get_text().encode('utf-8')

	start = t.find('Title')
	title_snp = ''
	num = start + len('Title') + 1
	letter = t[start + len('Title') + 1]
	while(letter != '\n'):
		title_snp = title_snp + letter
		num+=1
		letter = t[num]
	start = t.find('Characteristics')

	characteristic_snp = ''
	num = start + len('Characteristics\nprimary site: ')
	letter = t[start + len('Characteristics\nprimary site: ')]
	while(letter != '\n'):
		characteristic_snp = characteristic_snp + letter
		num+=1
		letter = t[num]
		if t[num:num+len('histology')] == 'histology':
			break

	temp = headers_gene[j][1:(len(headers_gene[j])-1)]
	link1 = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=" + temp
	page = urllib2.urlopen(link1)
	soup = BeautifulSoup(page)
	t =  soup.get_text().encode('utf-8')

	start = t.find('Title')
	title_gene = ''
	num = start + len('Title') + 1
	letter = t[start + len('Title') + 1]
	while(letter != '\n'):
		title_gene = title_gene + letter
		num+=1
		letter = t[num]
	#print title_snp
	start = t.find('Characteristics')

	characteristic_gene = ''
	num = start + len('Characteristics\nprimary site: ') 
	letter = t[start + len('Characteristics\nprimary site: ') ]
	while(letter != '\n'):
		characteristic_gene = characteristic_gene + letter
		num+=1
		letter = t[num]
		if t[num:num+len('histology')] == 'histology':
			break
	#print characteristic_gene


	if title_snp == title_gene and characteristic_snp == characteristic_snp:
		if (title_snp.replace('-','') +'_' + characteristic_snp.upper()) in lines:
			print "final: " + title_snp.replace('-','') +'_' + characteristic_snp.upper()
			index = lines.index(title_snp.replace('-','') +'_' + characteristic_snp.upper())
			final_cell_lines.append(title_snp.replace('-','') +'_' + characteristic_snp.upper())
			final_ics.append(ics[index])
			final_snp_data = np.append(final_snp_data, np.array(result_snp[i]).reshape(1,22419),axis=0)
			final_gene_data = np.append(final_gene_data, np.array(result_gene[j]).reshape(1,18926),axis=0)

		j+=1


np.save('temp_output/final_snp_data.npy', final_snp_data)
np.save('temp_output/final_gene_data.npy', final_gene_data)
np.save('temp_output/final_ics.npy', np.array(final_ics))
np.save('temp_output/final_cell_lines.npy', np.array(final_cell_lines))


