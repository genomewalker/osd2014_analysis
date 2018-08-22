#!/usr/bin/env python

import sys
import pandas as pd
from sklearn.cluster import AffinityPropagation

infile = sys.argv[1]
outfile = sys.argv[2]
damping = sys.argv[3]
maxit = sys.argv[4]
covit = sys.argv[5]
preference = sys.argv[6]

m = pd.read_csv(infile)
m_names = m.columns


if preference == 'None':
    af = AffinityPropagation(damping= float(damping), max_iter=int(maxit), convergence_iter=int(covit), copy=True, preference=None, affinity='precomputed').fit(m)
else:
    preference = float(preference)
    af = AffinityPropagation(damping= float(damping), max_iter=int(maxit), convergence_iter=int(covit), copy=True, preference=preference, affinity='precomputed').fit(m)





cluster_centers_indices = af.cluster_centers_indices_
labels = af.labels_

ap_f = open(outfile, 'w')
for x in range(len(labels)):
    ap_f.write(str(m_names[x])+"\t"+str(labels[x]) +"\n")
#        print(out, file=ap_f)
ap_f.close()

