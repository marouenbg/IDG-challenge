from Bio.ExPASy import ScanProsite
import numpy as np
import time

protStructure=np.load("../data/protStructure.npy").item()
kinaseDict={}
ATPDict   ={}

count=0
missing=0
for key, sequence in protStructure.items():
    count=count+1
    print(count)
    if count % 50 ==0:
        time.sleep(60)# sleep 1 mn for very 50 query to avoid timeout
    handle = ScanProsite.scan(seq=sequence)
    result = ScanProsite.read(handle)
    for i in range(len(result)): #I am looping over all results but there should be only one that has kinase
        if result[i]['signature_ac']=='PS50011':# Protein kinase domain
            kinaseDict[key]=sequence[result[i]['start']:result[i]['stop']]
        elif result[i]['signature_ac']=='PS00107':# ATP binding pocket
            ATPDict[key]=sequence[result[i]['start']:result[i]['stop']]
            
print(missing)
np.save('kinaseDomain.npy', kinaseDict)
np.save('ATPDomain.npy', ATPDict)
