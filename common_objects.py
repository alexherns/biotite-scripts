"""Creates a scaffold2bin dictionary object from the respective file"""
def scaf2bin_dictionary(file_handle, header=True, sep="\t"):
    if header==True:
        file_handle.readline()
    scaf2bin= {}
    for line in file_handle:
        scaffold, Bin= line.split(sep)[:2]
        scaf2bin[scaffold]= Bin
    return scaf2bin
