import pubchempy as pcp #api query
import pandas as pd #data loading
import os #to change dir

def isNaN(num):
    return num != num

#read drugcommons db
os.chdir("../data")
df=pd.read_csv("DtcDrugTargetInteractions.csv", usecols=['standard_inchi_key'],dtype={'standard_inchi_key': object})
print(df.shape)

#declare feature lists
smiles_list=[]

#fetch canonical smiles
type='inchikey' #query by inchi key

prev=[]
missing=0 #missing compounds
for key in df["standard_inchi_key"]:
	#skip repeated compounds
	if key==prev:
		continue
	prev=key 
	if not isNaN(key):
		results = pcp.get_compounds(key, type)
		if results==[]:
			#for missing compounds we can search using another key such as name
			missing=missing+1
			continue
		#take the first result here, should be fine because smiels are unique I think
		smiles_list.append(results[0].to_dict(properties=['canonical_smiles']))

print(smiles_list)

