import streamlit as st
# import importlib
import astromol
# astromol = importlib.import_module("astromol")
st.set_page_config(page_title="Astromol", layout="wide")

def main():
    st.write(f"version {astromol.__version__}")
    about_page()

def about_page():
    st.divider()
    with st.sidebar:
        """
            ## About _astromol_
            

            A Database of Molecules Detected in Space
            
            The main python library is created and mainted by Dr. Brett A. McGuire
            
            [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
            
            [Original source code](https://github.com/bmcguir2/astromol)
            
            [Read paper](https://iopscience.iop.org/article/10.3847/1538-4365/aae5d2)
            
            If you use `astromol` for your own work, please cite the Zenodo entry: [![DOI](https://zenodo.org/badge/360663606.svg)](https://zenodo.org/badge/latestdoi/360663606)
        """
        
if __name__ == "__main__":
    st.title("Astromol")
    # st.write("This web interface is created by A.N. Marimuthu.")
    main()
    