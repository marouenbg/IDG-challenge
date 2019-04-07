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
from Bio.ExPASy import ScanProsite

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
os.chdir("../data/Round2")
df=pd.read_csv("round_2_template.csv")
print(df.shape)

uniqueInchiKeys=df["Compound_InchiKeys"].unique() 

####2.Fetch chemical structures AND other features
#declare feature lists
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
for key in uniqueInchiKeys:
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
np.save("Round2_cactvsFingerprint.npy", cactvsFingerprint_dict)
#save key names: these are the inchi key of drugs (rows)
#the columns names are the featrues name (otehrFetaures variable)
with open('Round2_drug_names.txt', 'w') as f:
    for key in keys_list:
        f.write("%s\n" % key)
#save the other features matrix
chemFeatMat.to_csv("Round2_chemFeatMat.csv",sep=',')

####3.Fetch protein sequence
seq_dict={}
uniqueProts=df["UniProt_Id"].unique() 

missingProts=0
for prot in uniqueProts:
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

print(missingProts)
np.save('Round2_protStructure.npy', seq_dict)
####4.Extract kinase domain and ATP binding pocket
kinaseDict={}
ATPDict   ={}

amissing,kmissing,count=0,0,0
for key, sequence in seq_dict.items():
	count=count+1
	print(key)
	if count % 50 ==0:
		time.sleep(60) # sleep 1 mn for very 50 query to avoid timeout
	handle = ScanProsite.scan(seq=sequence)
	result = ScanProsite.read(handle)
	kinase,atp=0,0
	for i in range(len(result)): #I am looping over all results but there should be only one that ha$
		if result[i]['signature_ac']=='PS50011': # Protein kinase domain
			kinaseDict[key]=sequence[result[i]['start']:result[i]['stop']]
			kinase=1
		elif result[i]['signature_ac']=='PS00107': # ATP binding pocket
			ATPDict[key]=sequence[result[i]['start']:result[i]['stop']]
			atp=1
	if kinase==0:
		kmissing=kmissing+1
		print('kinase missing')
	if atp==0:
		amissing=amissing+1
		print('atp missing')

print(kmissing)
print(amissing)
np.save('Round2_kinaseDomain.npy', kinaseDict)
np.save('Round2_ATPDomain.npy', ATPDict)
