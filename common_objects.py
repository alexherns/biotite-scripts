import numpy as np

"""Creates a scaffold2bin dictionary object from the respective file"""
def scaf2bin_dictionary(file_handle, header=True, sep="\t"):
    if header==True:
        file_handle.readline()
    scaf2bin= {}
    for line in file_handle:
        scaffold, Bin= line.split(sep)[:2]
        scaf2bin[scaffold]= Bin
    return scaf2bin

"""Reads an esom.lrn file into a numpy array"""
def lrn2np(file_handle):
	return np.asarray([[float(i) for i in line.strip().split()[1:]] for line in file_handle if "%" not in line])
	
"""Writes an esom.lrn file from a numpy array given an output file handle"""
def writelrn(X, file_handle):
	num_rows, num_cols= np.shape(X)
	print num_rows
	print num_cols
	file_handle.write('%\t{0}\n%\t{1}\n% key \t{2}\n'.format(str(num_rows), str(num_cols), '\t'.join(["D"+str(d) for d in range(num_cols)])))
	for i, line in enumerate(X.tolist()):
		file_handle.write('{0}\t{1}\n'.format(str(i+1), '\t'.join([str(j) for j in line])))