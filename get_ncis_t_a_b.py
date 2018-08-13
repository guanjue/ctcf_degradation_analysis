import numpy as np
import matplotlib.pyplot as plt
import plotly.plotly as py

def get_ncis_t_a_b(matrix1, matrix2, output):
	data1 = open(matrix1,'r')
	data2 = open(matrix2,'r')

	### skip header
	data1.readline()
	data2.readline()

	d1_tri1 = {}
	d1_tri2 = {}
	for rec1,rec2 in zip(data1, data2): ### loop matrix at x axis
		d1 = [ x.strip() for x in rec1.split('\t') ]
		d2 = [ x.strip() for x in rec2.split('\t') ]
		### add 1 for each bin
		r1 = int(float(d1[1])*10) +1
		r2 = int(float(d2[1])*10) +1
		### get t value
		t = r1+r2

		### get t dict for sig1 & sig2
		if t in d1_tri1:
			d1_tri1[t].append( r1 )
			d1_tri2[t].append( r2 )
		else:
			d1_tri1[t] = [ r1 ]
			d1_tri2[t] = [ r2 ]
	data1.close()
	data2.close()

	print('read table DONE')
	### get t list
	d1_tri_sum = d1_tri1.keys() ### get t for NCIS
	t_list = np.unique(d1_tri_sum) ### get unique t 
	### sort t list
	t_order = np.argsort(t_list) ### sort t id
	t_list_sort = t_list[t_order] ### sort t

	### get r list
	d1_tri_tsum =[]
	for t in t_list_sort:
		### sum signals of bins with the same t
		d1_tmp = np.sum(d1_tri1[t]) ### get matrix1 sum, when matrix1+matrix2 = ti
		d2_tmp = np.sum(d1_tri2[t]) ### get matrix1 sum, when matrix1+matrix2 = ti
		d1_tri_tsum.append( [ t, d1_tmp, d2_tmp ] )
	d1_tri_tsum = np.array(d1_tri_tsum)
	print('ncis table DONE')
	#print(d1_tri_tsum)
	result = open(output,'w')
	for records in d1_tri_tsum:
		result.write(str(records[0])+'\t'+str(records[1])+'\t'+str(records[2])+'\n')
	result.close()

	fig, ax = plt.subplots()
	ax.scatter(np.log10(d1_tri_tsum[:,0].astype(float)), np.log10((d1_tri_tsum[:,1].astype(float) / d1_tri_tsum[:,2].astype(float))) )
	plt.savefig(output+".ncis_tr_plot.png")
	#plot_url = py.plot_mpl(fig, filename="ncis_tr_plot.png")

###########################################
# time python get_atac_t_ncis.py -i sample1.txt -j sample2.txt -o sample12.tr.txt
import getopt
import sys

def main(argv):
	try:
		opts, args = getopt.getopt(argv,"hi:j:o:")
	except getopt.GetoptError:
		print 'python extract_statepair_position.py -i sample1.tab -j sample2.tab -o output_file.txt'
		sys.exit(2)

	for opt,arg in opts:
		if opt=="-h":
			print 'python extract_statepair_position.py -i sample1.tab -j sample2.tab -o output_file.txt'
			sys.exit()
		elif opt=="-i":
			matrix1=str(arg.strip())
		elif opt=="-j":
			matrix2=str(arg.strip())
		elif opt=="-o":
			output=str(arg.strip())
	get_ncis_t_a_b(matrix1, matrix2, output)

if __name__=="__main__":
	main(sys.argv[1:])

