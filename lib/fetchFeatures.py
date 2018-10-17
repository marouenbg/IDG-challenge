import pubchempy as pcp #api query
import pandas as pd #data loading
import os #to change dir
import requests #query uniprot
from Bio import SeqIO #parse fasta
from requests.adapters import HTTPAdapter #for retries
from requests.packages.urllib3.util.retry import Retry
import numpy as np #save dict
from scipy import sparse #for spare matrices
import math #to compute log10

def isNaN(num):
	return num != num

def requests_retry_session(
	###########
	#from: https://www.peterbe.com/plog/best-practice-with-retries-with-requests
	###########
	retries=3,
	backoff_factor=0.3,
	status_forcelist=(500, 502, 504),
    	session=None,
	):
	session = session or requests.Session()
	retry = Retry(
		total=retries,
		read=retries,
		connect=retries,
		backoff_factor=backoff_factor,
		status_forcelist=status_forcelist,
    		)
	adapter = HTTPAdapter(max_retries=retry)
	session.mount('http://', adapter)
	session.mount('https://', adapter)
	return session


#read drugcommons db
os.chdir("../data")
df=pd.read_csv("DtcDrugTargetInteractions.csv", usecols=['standard_inchi_key', 'target_id', 'standard_type', 'standard_relation', 'standard_value', 'standard_units'],\
	dtype={'standard_inchi_key': object, 'target_id':object})
print(df.shape)

####1.Fetch pkd values
#create pkd matrix that has drugs as colmuns and proteins as rows
#Unique compounds
uniqueDrugs=df["standard_inchi_key"].unique()
uniqueDrugs=uniqueDrugs[1:]
#Unique protein IDs
#start from 1 because first entry is nan
uniqueProt=df["target_id"].unique() #take the first code for now (6 first letters)
uniqueProt = [prot[:6] for prot in uniqueProt[1:]] #take out na which is the first occurence

colnames=uniqueDrugs
rownames=uniqueProt

rowId,colId,data=[],[],[]
for i in range(df.shape[0]):
        prot=df["target_id"].iloc[i]
        drug=df["standard_inchi_key"].iloc[i]
        #only take pkd assays in NM because in the webinar, they said that assays reported in % are not reliable, will get back to that
        if not isNaN(drug) and not isNaN(prot) and df["standard_type"].iloc[i]=='KDAPP' and df["standard_units"].iloc[i]=="NM":
                prot=prot[:6] # take the first protein
                rowId.append( rownames.index(prot) )
                colId.append( np.where(colnames==drug)[0][0] )
                # I did not consider values that repeat; maybe average them
                #assume equality for now
                data.append(-math.log(df["standard_value"].iloc[i],10))

#build sparse matrix (0 can be treated as a missing value)
pkdMat = sparse.coo_matrix((data, (rowId, colId)))
#save (4951 interaction rows 317, cols 213)
sparse.save_npz('pkdMat.npz', pkdMat)
print(len(rowId))
print(len(colId))

####2.Fetch chemical structures
#declare feature lists
smiles_dict={}

#fetch canonical smiles
type='inchikey' #query by inchi key

missingDrugs=0 #missing compounds
for key in uniqueDrugs[np.unique(colId)]:
	print(key)
	if not isNaN(key):
		results = pcp.get_compounds(key, type)
		if results==[]:
			#for missing compounds we can search using another key such as name
			missingDrugs=missingDrugs+1
			print('missing')
			continue
		#take the first result here, should be fine because smiles are unique I think
		smiles=results[0].to_dict(properties=['canonical_smiles'])
		smiles_dict.update({key:smiles['canonical_smiles']})


#print(smiles_dict)
print(missingDrugs)
np.save('chemStructure.npy', smiles_dict) 

####2.Fetch protein sequence
seq_dict={}

missingProts=0
for prot in np.array(uniqueProt)[np.unique(rowId)]:
	if not isNaN(prot):
		print(prot)
		r=requests_retry_session().get('https://www.uniprot.org/uniprot/' + prot + '.fasta', timeout=10)
		if r.status_code != requests.codes.ok:
			missingProts=missingProts+1
			print('missing')
			continue
		#parsing fasta with biopython requeries a file
		text_file = open("prot.fasta", "w")
		text_file.write(r.text)
		text_file.close()

		for record in SeqIO.parse("prot.fasta", "fasta"):
			seq_dict.update({prot:record.seq})
		#deleting file
		os.remove("prot.fasta")

print(missingProts) #1 records are missing
np.save('protStructure.npy', seq_dict) 

