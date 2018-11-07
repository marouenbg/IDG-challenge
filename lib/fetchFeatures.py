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
import time #to sleep

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

rowId,colId,data,relData=[],[],[],[]
for i in range(df.shape[0]):
	prot=df["target_id"].iloc[i]
	drug=df["standard_inchi_key"].iloc[i]
	#only take pkd assays in NM because in the webinar, they said that assays reported in % are not reliable, will get back to that
	if not isNaN(drug) and not isNaN(prot) and df["standard_type"].iloc[i] in ['Kd','KD','KDAPP'] and df["standard_units"].iloc[i]=="NM":
		if(df["standard_value"].iloc[i]) == 0:
			continue #to avoid log error // maybe save in kd not pkd 
		prot=prot[:6] # take the first protein
		rowId.append( rownames.index(prot) )
		colId.append( np.where(colnames==drug)[0][0] )
		# I did not consider values that repeat; maybe average them
		#assume equality for now
		data.append(-math.log(df["standard_value"].iloc[i],10))
		if df["standard_relation"].iloc[i] == '=':
			rel=1
		elif df["standard_relation"].iloc[i] == '<':
			rel=2
		else:
			rel=3
		relData.append(rel)

#build sparse matrix (0 can be treated as a missing value)
pkdMat = sparse.coo_matrix((data, (rowId, colId)))
relMat = sparse.coo_matrix((relData, (rowId, colId)))
#save (4951 interaction rows 317, cols 213)
sparse.save_npz('pkdMat.npz', pkdMat)
sparse.save_npz('relMat.npz', relMat)
print(len(rowId))
print(len(colId))

####2.Fetch chemical structures AND other features
#declare feature lists
#these are the drug smiles dictioanry using the inchi keys as keys
smiles_dict={}
#this is the cactvs fingerprint dictionary using the inchi keys as keys
cactvsFingerprint_dict={}

#fetch canonical smiles
type='inchikey' #query by inchi key
otherFeaturesName=["atom_stereo_count", "bond_stereo_count", "charge","complexity", \
"covalent_unit_count","defined_atom_stereo_count", "exact_mass", "h_bond_acceptor_count", \
"h_bond_donor_count","heavy_atom_count","isotope_atom_count", "molecular_weight", \
"monoisotopic_mass","rotatable_bond_count", "tpsa", "xlogp"]

#Int: bonds, atoms, elements, fingerprint, molecular formula, record

#the following dataframe contains a set of nuelrical features collected from pubchem
#the columns are features from the otehrFeaturesName list
#the rowas are the drugs listed in keys_list variable by inchi keys
chemFeatMat= pd.DataFrame(columns=otherFeaturesName)
missingDrugs=0 #missing compounds
counter=0
keys_list=list()
for key in uniqueDrugs[np.unique(colId)]:
	print(key)
	if not isNaN(key):
		try: 
			results = pcp.get_compounds(key, type)
		except Exception:
			time.sleep(60)
			results = pcp.get_compounds(key, type) # try again
		if results==[]:
			#for missing compounds we can search using another key such as name
			missingDrugs=missingDrugs+1
			print('missing')
			continue
		counter=counter+1
		#take the first result here, should be fine because smiles are unique I think
		smiles=results[0].to_dict(properties=['canonical_smiles'])
		smiles_dict.update({key:smiles['canonical_smiles']})
		#cactvs fingerprint
		cactvs=results[0].to_dict(properties=["fingerprint"])
		cactvsFingerprint_dict.update({key:cactvs['fingerprint']})
		#other features
		otherFeat=results[0].to_dict(properties=otherFeaturesName)
		chemFeatMat.loc[counter]=np.zeros(len(otherFeaturesName))
		chemFeatMat.loc[counter]=otherFeat.values()
		#keys name
		keys_list.append(key)

#print(smiles_dict)
print(missingDrugs) #414
np.save('chemStructure.npy', smiles_dict) 
np.save("cactvsFingerprint.npy", cactvsFingerprint_dict)
#save key names: these are the inchi key of drugs (rows)
#the columns names are the featrues name (otehrFetaures variable)
with open('drug_names.txt', 'w') as f:
    for key in keys_list:
        f.write("%s\n" % key)
#save the other features matrix
chemFeatMat.to_csv("chemFeatMat.csv",sep=',')

####3.Fetch protein sequence
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


