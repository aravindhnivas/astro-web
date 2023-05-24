import streamlit as st
import UCLCHEM as uclchem
from pages.UCLCHEM.parameters import (
    get_parameters, get_behavioural_parameters, 
    get_input_output_parameters, get_integration_controls
)

def about_page():
    st.divider()
    with st.sidebar:
        """
            ## About _UCLCHEM_

            A Gas-Grain Chemical Code for astrochemical modelling
            
            UCLCHEM is a gas-grain chemical code for astrochemical modelling that can be used as a stand alone Fortran program or a Python module. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas.
            
            visit [here](https://uclchem.github.io/) for more information.
            
            [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
            
            Code and models are available on [GitHub](https://github.com/uclchem/UCLCHEM), 
            and  [published paper](https://ui.adsabs.harvard.edu/abs/2017AJ....154...38H/abstract) 📃
            
        """

def simple_model():
    st.subheader("Simple model")
    st.write("The simple model is a 1D model of a static cloud with a constant density and temperature. The model is run for a specified time and the abundances of the chemical species are outputted at the end of the run.")


def main():
    
    st.header("UCLCHEM")
    st.write("A Gas-grain Chemical Code for Clouds, Cores, and C-Shocks")
    st.divider()
    
    st.subheader("Fine-tuning the model")
    
    with st.expander("Parameters"):
        parameters = get_parameters()
    with st.expander("Behavioural parameters"):
        behaviour_parameters = get_behavioural_parameters()
    with st.expander("Input and Output parameters"):
        input_output_parameters = get_input_output_parameters()
    with st.expander("Integration Controls"):
        integration_controls = get_integration_controls()
    
    # st.write(parameters, behaviour_parameters, input_output_parameters, integration_controls)
    
    st.divider()
    
    st.markdown("""
        ## A Simple Cloud

        UCLCHEM's `cloud` model is a spherical cloud of isothermal gas. We can keep a constant density or have it increase over time following a freefall equation. This model is generally useful whenever you want to model a homogeneous cloud of gas under constant conditions. For example, in the inner parts of a molecular cloud where Av $\gtrsim$ 10 there are very few depth dependent processes. You may wish to model the whole of this UV shielded portion of the cloud with a single `cloud` model.

        Due to the large number of parameters in a chemical model and the way fortran and python interaction, we find it is easiest to do parameter input through python dictionaries. In this block, we define param_dict which contains the parameters we wish to modify for this run. Every `uclchem.model` function accepts a dictionary as an optional argument. Every parameter has a default value which is overriden if that parameter is specified in this dictionary. You can find a complete list of modifiable parameters and their default values as above.
        
    """)
    
    species = st.text_input("outSpecies", value='SO, CO')
    out_species = [_.strip() for _ in species.split(',')]
    
    param_dict = parameters | behaviour_parameters | input_output_parameters | integration_controls
    # result = uclchem.model.cloud(param_dict=param_dict, out_species=out_species)
    
    st.write(uclchem)
    
if __name__ == "__main__":
    main()
    about_page()
        