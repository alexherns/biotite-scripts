import os
import socket

DATA3DIR= '/data3/borehole'

if socket.gethostname() == 'biotite':
    FULLSCAF2BINPATH= '/data3/borehole/scaffolds2bins/horonobe.'
    NRSCAF2BINPATH= '/data3/borehole/scaffolds2bins/horonobe_nr.scaf2bin.kch'
else:
    BIN_TO_BIN_CSV='/Users/alexh/Documents/berkeley/banfield_lab/organism_info/best_choices-LF.csv'
    PUBLICATION_NAMES='/Users/alexh/Documents/berkeley/banfield_lab/organism_info/publication_names.csv'
    FULLSCAF2BINPATH=\
            '/Users/alexh/Documents/berkeley/banfield_lab/scaffolds2bins/horonobe.scaffolds_to_bin.kch'
    NRSCAF2BINPATH= \
            '/Users/alexh/Documents/berkeley/banfield_lab/scaffolds2bins/horonobe_nr.scaf2bin.kch'
