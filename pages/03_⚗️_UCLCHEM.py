import streamlit as st
from UCLCHEM.src import uclchem
from pathlib import Path as pt
import re
from pages.UCLCHEM.parameters import (
    get_parameters, get_behavioural_parameters, 
    get_input_output_parameters, get_integration_controls
)
from time import perf_counter
import pandas as pd
st.set_page_config(layout='wide')
loc = pt("./pages/UCLCHEM/outputs").absolute()


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

def set_loc(filename: str):
    if not filename.endswith('.dat'):
        filename = filename + ".dat"
    
    if not loc.exists(): 
        loc.mkdir()
        
    return str(loc / filename)

    
def main():
    
    st.header("UCLCHEM")
    st.write("A Gas-grain Chemical Code for Clouds, Cores, and C-Shocks")
    # st.write(f"version: {uclchem.__version__}")
    st.divider()
    
    st.subheader("Fine-tuning the model")
    
    with st.expander("Parameters"):
        parameters = get_parameters()
    with st.expander("Behavioural parameters"):
        behaviour_parameters = get_behavioural_parameters()
    with st.expander("Input and Output parameters"):
        input_output_parameters = get_input_output_parameters()
        
        input_output_parameters_filtered = {
            key: set_loc(value)
            for key, value in input_output_parameters.items() if key != 'writeStep' and value
        }
        input_output_parameters_filtered['writeStep'] = input_output_parameters['writeStep']
        
        
    with st.expander("Integration Controls"):
        integration_controls = get_integration_controls()
    
    st.divider()
    
    st.markdown("""
        ## A Simple Cloud

        UCLCHEM's `cloud` model is a spherical cloud of isothermal gas. We can keep a constant density or have it increase over time following a freefall equation. This model is generally useful whenever you want to model a homogeneous cloud of gas under constant conditions. For example, in the inner parts of a molecular cloud where Av $\gtrsim$ 10 there are very few depth dependent processes. You may wish to model the whole of this UV shielded portion of the cloud with a single `cloud` model.

        Due to the large number of parameters in a chemical model and the way fortran and python interaction, we find it is easiest to do parameter input through python dictionaries. In this block, we define param_dict which contains the parameters we wish to modify for this run. Every `uclchem.model` function accepts a dictionary as an optional argument. Every parameter has a default value which is overriden if that parameter is specified in this dictionary. You can find a complete list of modifiable parameters and their default values as above.
        
    """)
    
    species = st.text_input("outSpecies", value='SO, CO+')
    
    out_species = [_.strip() for _ in species.split(',')]
    
    param_dict = parameters | input_output_parameters_filtered | behaviour_parameters | integration_controls
    time_start = perf_counter()
    
    if st.button('Run calculations'):
        
        for mol in out_species:
            if not re.match(r'^[a-zA-Z][a-zA-Z0-9+,\-]+$', mol):
                st.error(f"{mol} contains invalid characters. Please change it to continue")
                return
            
        result = uclchem.model.cloud(param_dict=param_dict, out_species=out_species)
        # st.write(result)
        
        status = result[0]
        if status > 0:
            st.success(f'Finished in {(perf_counter() - time_start):.2f} seconds')
            st.markdown("### Final abundances")
            
            final_abundances_of_species = pd.DataFrame({
                "name": [mol for mol in out_species],
                "abundance": [f"{abundance:.2e}" for abundance in result[1:]]
            })
            st.dataframe(final_abundances_of_species, use_container_width=True)
            
            result_df: pd.DataFrame = None
            file_saved = 'outputFile' in input_output_parameters_filtered and input_output_parameters_filtered['outputFile']
            
            with st.expander("Show full output"):
                if file_saved:
                    result_df = uclchem.analysis.read_output_file(input_output_parameters_filtered['outputFile'])
                    st.dataframe(result_df)
                else:
                    st.warning("No outputFile mentioned so cannot continue with further analysis")
            
            if file_saved:
                   
                with st.expander("Check: elemental conservation", expanded=True):
                    # st.markdown("### Elemental conservation")
                    st.markdown("""
                        To test whether we conserve elements. Each entry gives the change in the total abundance of an element as a percentage of the original abundance. In an ideal case, these values are 0\% indicating the total abundance at the end of the model is exactly the same as the total at the start.
                        
                        Changes of less than 1\% are fine for many cases but if they are too high, you could consider changing the `reltol` and `abstol` parameters that control the integrator accuracy. They are error tolerance so smaller values lead to smaller errors and (usually) longer integration times. The default values were chosen by running a large grid of models and choosing the tolerances with the lowest average run time from those that conserved elements well and rarely failed. Despite this, there are no one-size-fits-all perfect tolerances and you may run into issues with different networks or models.
                                    
                    """)
                    
                    element_list = st.text_input('element_list', value='H, N, C, O, S')
                    element_list_arr = [_.strip() for _ in element_list.split(',')]
                    conservation=uclchem.analysis.check_element_conservation(result_df, element_list=element_list_arr)
                    st.write("Percentage change in total abundances:")
                    st.dataframe(conservation)
                
                with st.expander("Plotting result", expanded=True):
                    st.markdown("""
                        Note the use of $ symbols in the species list below, this gets the total ice abundance of a species. For two phase models, this is just the surface abudance but for three phase it is the sum of surface and bulk.
                    """)
                    
                    species_list = st.text_input('species_list', value='H, H2, $H, $H2, H2O, $H2O, CO, $CO, $CH3OH, CH3OH')
                    species_list_arr = [_.strip() for _ in species_list.split(',')]
                    
                    col1, col2, col3, col4 = st.columns(4)
                    figwidth = col1.number_input('width', value=10)
                    figHeight = col2.number_input('height', value=7)
                    
                    xlim = col3.text_input('xlim', value='1e3, 1e6')
                    ylim = col4.text_input('ylim', value='1e-15, 1')
                    
                    fig, ax = uclchem.analysis.create_abundance_plot(result_df, species_list_arr, figsize=(figwidth, figHeight))
                    
                    ax = ax.set(
                        xscale="log", 
                        ylim=[float(_) for _ in ylim.split(',')], 
                        xlim=[float(_) for _ in xlim.split(',')]
                    )
                    
                    st.pyplot(fig)
                
        else:
            st.error('Error occured during calculation')
            
if __name__ == "__main__":
    
    main()
    about_page()
    