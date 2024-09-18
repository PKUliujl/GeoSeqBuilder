
import pandas as pd
import numpy as np

import random
import os
dataR = pd.read_csv(os.path.join(os.path.dirname(os.path.abspath(__file__)),'AA_rotamer_list'),header=None,sep='\s+')
interval = np.linspace(-180,180,49)
def samplingr(res,chi,chi_class):
    #dataR = pd.read_csv('/home/lhlai_pkuhpc/lustre3/liujl/RD/distribution/AA_rotamer_list',header=None,sep='\s+')
    #interval = np.linspace(-180,180,49)
    rotamer_list = dataR[ (dataR.iloc[:,0] == res) & (dataR.iloc[:,1] == chi) & (dataR.iloc[:,2] == chi_class) ].iloc[0,:].to_list()
    rn = random.choices([i for i in range(15)],weights = rotamer_list[3:], k=1)[0]
    mini_interval = np.linspace( interval[chi_class], interval[ chi_class+1], 16)
    angle = np.random.uniform(mini_interval[rn], mini_interval[rn+1],1)[0]
    return angle



'''
AA = "ARG LYS MET GLU GLN ASP ASN ILE LEU HIS TRP TYR PHE PRO THR VAL SER CYS"

for res in AA.split():
    print(sampleingir(res,0,47))
'''
