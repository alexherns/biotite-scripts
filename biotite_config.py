import os
import socket

DATA3DIR= '/data3/borehole'

if socket.gethostname() == 'biotite':
    SCAF2BINPATH= '/data3/borehole/scaffolds2bins/horonobe_nr.scaf2bin.kch'
else:
    SCAF2BINPATH= '/Users/alexh/Documents/berkeley/banfield_lab/scaffolds2bins/horonobe_nr.scaf2bin.kch'
