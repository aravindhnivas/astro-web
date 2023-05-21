import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw

from joblib import load
import numpy as np
from umda.data import load_pipeline

import pandas as pd

st.title("Unsupervised Molecule Discovery in Astrophysics")
st.write("This is a web app for the UMDA project.")
st.write('Read more on : Lee+, ‘Machine Learning of Interstellar Chemical Inventories’, ApJL, vol. 917, no. 1, p. L6, Aug. 2021, doi: 10.3847/2041-8213/ac194b.')

st.header('TMC-1 abundances prediction')



# load a wrapper class for generating embeddings
embedder = load_pipeline()
regressors = load("models/regressors.pkl")
regressor_list = list(regressors.keys())



smiles = st.text_input('Enter SMILES string', 'C1=C=C=C=C=C1, C=O')
mol_smiles_lists = [mol.strip() for mol in smiles.split(",")]
# print(smiles.split(","))

chem_images = []
for mol in mol_smiles_lists:
    molecule = Chem.MolFromSmiles(mol.strip())
    img = Draw.MolToImage(molecule)
    chem_images.append(img)

st.image(chem_images, caption=mol_smiles_lists)

vecs = np.vstack([embedder.vectorize(smi) for smi in mol_smiles_lists])

regressor_model_lists = st.multiselect("Choose a regressor model", regressor_list, default=["gbr", "svr", "rfr"])

predictions = []
for regressor_model_name in regressor_model_lists:
    regressor = regressors.get(regressor_model_name)
    prediction = regressor.predict(vecs)
    predictions.append(prediction)
    
# print(predictions)

predictions = np.vstack(predictions).T
# st.write(predictions)
# print(predictions)

# st.subheader(f"Predicted abundances")

df = pd.DataFrame(predictions, columns=regressor_model_lists, index=mol_smiles_lists)
st.table(df)
