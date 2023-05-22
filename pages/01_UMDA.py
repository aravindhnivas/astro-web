import streamlit as st

from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors

from joblib import load
import numpy as np
# import importlib
import pandas as pd
# import streamlit.components.v1 as components
from umda.data import load_pipeline
# load_pipeline = importlib.import_module("umda.data").load_pipeline
st.set_page_config(page_title="UMDA", layout="wide")


@st.cache_data
def get_regressor():
    # load a wrapper class for generating embeddings
    embedder = load_pipeline()
    regressors = load("models/regressors.pkl")
    regressor_list = list(regressors.keys())
    return embedder, regressors, regressor_list

def main():
    st.header('Column density prediction towards TMC-1')
    st.warning("Enter a SMILES string of a molecule in the sidebar and choose a regressor model(s) to predict its column density. Some examples are provided below.")
    
    st.divider()
    
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
        
    st.subheader("Chemical structures of the molecules")
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

    st.subheader("Predicted Column Density (log$_{10}$ cm$^{-2}$) for the selected models")
    
    df = pd.DataFrame(predictions, columns=regressor_model_lists, index=formula_lists)
    st.dataframe(df, use_container_width=True)
    
    about_page()

def about_page():
    
    with st.sidebar:
        st.divider()
        """
            ## About _UMDA_
            
            The UMDA project is created by [Dr. Kelvin Lee](https://laserkelvin.github.io/).
            
            Code and models are available on [GitHub](https://github.com/laserkelvin/umda), and published in [Ap. J. Letters](https://iopscience.iop.org/article/10.3847/2041-8213/ac194b)
            
            If you used the list of recommendations generated from this work as part of your own observations or work, please cite the Zenodo entry: [![DOI](https://zenodo.org/badge/360663606.svg)](https://zenodo.org/record/5080543)
        """
    
if __name__ == "__main__":
    st.title("Unsupervised Molecule Discovery in Astrophysics (UMDA)")
    # st.write("This web interface is created by A.N. Marimuthu.")
    st.divider()
    
    main()
