import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw, MolToSmiles, rdMolDescriptors

from joblib import load
import numpy as np
import importlib
import pandas as pd
# import streamlit.components.v1 as components
load_pipeline = importlib.import_module("umda.data").load_pipeline
st.set_page_config(layout="wide")

@st.cache_data
def get_regressor():
    # load a wrapper class for generating embeddings
    embedder = load_pipeline()
    regressors = load("models/regressors.pkl")
    regressor_list = list(regressors.keys())
    return embedder, regressors, regressor_list

def main():
    # Create a slider widget in the side panel
    # with st.sidebar:
    #     slider_val = st.slider('Select a value', 0, 100, 50)
    
    st.title("Unsupervised Molecule Discovery in Astrophysics (UMDA)")
    st.write("This is a web interface for the UMDA project created by A.N. Marimuthu.")
    st.write('Read more on : Lee+, ‘Machine Learning of Interstellar Chemical Inventories’, ApJL, vol. 917, no. 1, p. L6, Aug. 2021, doi: 10.3847/2041-8213/ac194b.')

    st.header('Column density prediction towards TMC-1')
    
    with st.sidebar:
        # mol_representation = st.radio("Choose a molecular representation", ("SMILES", "Chemical formula"))
        mol_representation = "SMILES"
        molecules_name = st.text_input(f'Enter molecule in "{mol_representation}" format', 'CCC#N, C#CC#[O+], CC1CCCCC1, Cc1ccccc1')
        if molecules_name == "":
            st.warning("Please enter SMILES string of a molecule to predict its column density.")
            return 
    
    molecules_name_lists = [mol.strip() for mol in molecules_name.split(",")]
    formula_lists = []
    chem_images = []
    smiles_lists = []
    
    for mol in molecules_name_lists:
        
        smiles = mol
        # Generate an RDKit molecule object from the SMILES string
        mol_obj = Chem.MolFromSmiles(smiles)
        # Generate a chemical formula from the molecule object
        formula = rdMolDescriptors.CalcMolFormula(mol_obj)
        
        formula_lists.append(formula)
        smiles_lists.append(smiles)
        
        img = Draw.MolToImage(mol_obj)
        chem_images.append(img)
        
    st.image(chem_images, caption=formula_lists)

    embedder, regressors, regressor_list = get_regressor()
    vecs = np.vstack([embedder.vectorize(smi) for smi in smiles_lists])

    with st.sidebar:
        regressor_model_lists = st.multiselect("Choose a regressor model", regressor_list, default=["gbr", "svr", "rfr"])

    if regressor_model_lists == []:
        st.warning("Please select at least one model.")
        return
    
    predictions = []
    for regressor_model_name in regressor_model_lists:
        regressor = regressors.get(regressor_model_name)
        prediction = regressor.predict(vecs)
        predictions.append(prediction)
        

    predictions = np.vstack(predictions).T

    st.subheader("Predicted Column Density (log$_{10}$ cm$^{-2}$) towards TMC-1 using the selected models")
    df = pd.DataFrame(predictions, columns=regressor_model_lists, index=formula_lists)
    st.table(df)
    
if __name__ == "__main__":
    main()
