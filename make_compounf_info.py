
import pubchempy as pcp
from rdkit.Chem import AllChem
from rdkit import Chem
import numpy as np
import pickle

# load file
train_set = []
with open('data/TrainSet.txt') as file:
    for f in file:
        line = f.split('\t')
        line[-1] = line[-1].split('\n')[0]
        train_set.append(line)
        
my_columns = train_set[0]

import pandas as pd

df = pd.DataFrame(train_set[1:],columns = my_columns)
id_list = df['Compound Identifier'].unique().tolist()

properties = ['IUPACName', 'MolecularFormula', 'MolecularWeight', 'XLogP', 'TPSA', 'CanonicalSMILES']

#get informations from pubChem
chem_infos = []
for cid in id_list:
    chem_info = pcp.get_properties(properties,cid)
    chem_infos.append(chem_info)

#ECFP4
radius = 2
nBits = 4096
morgan_fp = []
for info_list in chem_infos:
    info = info_list[0]
    mol = Chem.MolFromSmiles(info['CanonicalSMILES'])
    fp = [i for i in AllChem.GetMorganFingerprintAsBitVect(mol, radius, nBits)]
    morgan_fp.append(fp)
morgan_fp = np.array(morgan_fp)

compound_dic = {}
i = 0
for info_list in chem_infos:
    info = info_list[0]
    compound_dic[info['CID']] = morgan_fp[i]
    i += 1

#ファイル書き込み
f = open('data.45/compound_dic.pickle',mode='wb')
pickle.dump(compound_dic,f)
f.close()