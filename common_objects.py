import numpy as np

"""Function: scaf2bin_dictionary -- Creates a scaffold2bin dictionary object from the respective file"""
def scaf2bin_dictionary(file_handle, header=True, sep="\t"):
    if header==True:
        file_handle.readline()
    scaf2bin= {}
    for line in file_handle:
        scaffold, Bin= line.split(sep)[:2]
        scaf2bin[scaffold]= Bin
    return scaf2bin

"""Function: lrn2np -- Reads an esom.lrn file into a numpy array"""
def lrn2np(file_handle):
	return np.asarray([[float(i) for i in line.strip().split()[1:]] for line in file_handle if "%" not in line])
	
"""Function: writelrn -- Writes an esom.lrn file from a numpy array given an output file handle"""
def writelrn(X, file_handle):
	num_rows, num_cols= np.shape(X)
	print num_rows
	print num_cols
	file_handle.write('%\t{0}\n%\t{1}\n% key \t{2}\n'.format(str(num_rows), str(num_cols), '\t'.join(["D"+str(d) for d in range(num_cols)])))
	for i, line in enumerate(X.tolist()):
		file_handle.write('{0}\t{1}\n'.format(str(i+1), '\t'.join([str(j) for j in line])))

"""Function: motif2regex -- Converts a protein motif query to a regex search"""
def motif2regex(motif):
	return motif.replace('X', '.').replace('x', '.').replace('(', '{').replace(')', '}')

class UClustSeq():
    def __init__(self, line):
        self.cluster= int(line[1])
        self.query= line[8]

def uclust2dict(file_handle):
    from collections import defaultdict
    d= defaultdict(list)
    for line in file_handle:
        line= line.strip().split('\t')
        d[line[1]].append(line[8])
    return d
