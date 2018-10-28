# IDG-challenge
This is the code base for the submission of Let_Data_Talk team of the
[IDG-Dream kinase prediction challenge](https://www.synapse.org/#!Synapse:syn15667962/wiki/583305)

You need to install biopython

`pip3 install biopython`

Features:

- ATPDomain.npy     : dictionary of the sequence of ATP binding pockets by Uniprot ID

- kinaseDomain.npy  : dictionary of the sequence of kinase domain by Uniprot ID

- protStructure.npy : dictionary of the sequence of the whole protein by Uniprot ID

Dictionaries in numpy format are read as follows:

`dict=numpy.load("ATPDomain.npw").item()`

- chemStructure.npy : dictionary of the canonical smiles of the drugs by inchi keys

- pkdMat.npz        : sparse matrix of the pkd values between proteins (rows) drugs (columns)

- relMat.npz        : dataframe detailing the relation between the pkd and the value in pkdMat
		      1 means equal, 2 means less than, 3 means larger than
		      e.g., pkd[1,1]=4 and relmat[1,1]=2 means that the pkd value between drug 1
		      and protein 1 is less than 4.

		      The colnames of pkdMat and relMat are unique drug entries of the DTC databse
			`df["standard_inchi_key"].unique()`

		      The rownames of the pkdMat and relMat are unique protein entries of the DTC database	
			`df["target_id"].unique()`
