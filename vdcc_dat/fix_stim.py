#!/usr/bin/env python
from numpy import *

raw_dir = './vdcc_pq_data/'
files = ['VDCC_PQ_C01.dat',
         'VDCC_PQ_C10.dat',
         'VDCC_PQ_C12.dat',
         'VDCC_PQ_C21.dat',
         'VDCC_PQ_C23.dat',
         'VDCC_PQ_C32.dat',
         'VDCC_PQ_C34.dat',
         'VDCC_PQ_C43.dat',
         'VDCC_PQ_Ca.dat'
        ]

for fn in files:
   data = array([[float(x) for x in s.split()] for s in open(raw_dir + fn,'r').read().split('\n') if s != ''])
   tdata = [-data[2,0], 0] + data[2:93]
   tdata[len(tdata)-1,0] = 0.005
   outf = open(fn,'w')
   for t in tdata:
     outf.write('%g %g\n' % (t[0],t[1]))
   outf.close()
