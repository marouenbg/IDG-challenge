import pubchempy as pcp #api query
import pandas as pd #data loading
import os #to change dir
import requests #query uniprot
from Bio import SeqIO #parse fasta
from requests.adapters import HTTPAdapter #for retries
from requests.packages.urllib3.util.retry import Retry

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
df=pd.read_csv("DtcDrugTargetInteractions.csv", usecols=['standard_inchi_key', 'target_id'],dtype={'standard_inchi_key': object})
print(df.shape)

####1.Fetch chemical structures
#declare feature lists
smiles_list=[]

#fetch canonical smiles
type='inchikey' #query by inchi key

#Unique compounds
uniqueDrugs=df["standard_inchi_key"].unique()

missingDrugs=0 #missing compounds
for key in uniqueDrugs:
	if not isNaN(key):
		results = pcp.get_compounds(key, type)
		if results==[]:
			#for missing compounds we can search using another key such as name
			missingDrugs=missingDrugs+1
			continue
		#take the first result here, should be fine because smiels are unique I think
		smiles_list.append(results[0].to_dict(properties=['canonical_smiles']))


print(smiles_list)
print(missingDrugs)

####2.Fetch protein sequence
seq_list=[]

#Unique protein IDs
uniqueProt=df["target_id"].unique()
print(uniqueProt)
print(len(uniqueProt))

missingProts=0
for prot in uniqueProt:
	if not isNaN(prot):
		r=requests_retry_session().get('https://www.uniprot.org/uniprot/' + prot + '.fasta', timeout=100)
		if r.status_code != requests.codes.ok:
			missingProts=missingProts+1
			continue
		#parsing fasta with biopython requeries a file
		text_file = open("prot.fasta", "w")
		text_file.write(r.text)
		text_file.close()

		for record in SeqIO.parse("prot.fasta", "fasta"):
			seq_list.append(record.seq)
		#deleting file
		os.remove("prot.fasta")

print(seq_list)
print(missingProts)
