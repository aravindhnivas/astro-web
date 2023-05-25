import streamlit as st
from UCLCHEM.src import uclchem
from pathlib import Path as pt
# import re
from pages.UCLCHEM.parameters import (
    get_parameters, get_behavioural_parameters, 
    get_input_output_parameters, get_integration_controls
)

from time import perf_counter
import pandas as pd
pd.options.plotting.backend = "plotly"

# import plotly.express as px

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


def simple_cloud_model_calc():
    
    # if status < 0: return
    
    # result_df: pd.DataFrame = None
    # with st.expander("Show full output"):
    #     if file_saved:
    #         result_df = uclchem.analysis.read_output_file(param_dict['outputFile'])
    #         st.dataframe(result_df)
    #     else:
    #         st.warning("No outputFile mentioned so cannot continue with further analysis")
    
    if not pt(param_dict['outputFile']).exists():
        return st.warning("After setting appropriate parameters (above); Run calculation to show results")
    
    tab1, tab2 = st.tabs(['Results', 'Elemental conservation'])
    result_df = uclchem.analysis.read_output_file(param_dict['outputFile'])
    
    # with st.expander("Check: elemental conservation"):
    with tab2:
        st.markdown("""
            To test whether we conserve elements. Each entry gives the change in the total abundance of an element as a percentage of the original abundance. In an ideal case, these values are 0\% indicating the total abundance at the end of the model is exactly the same as the total at the start.
            
            Changes of less than 1\% are fine for many cases but if they are too high, you could consider changing the `reltol` and `abstol` parameters that control the integrator accuracy. They are error tolerance so smaller values lead to smaller errors and (usually) longer integration times. The default values were chosen by running a large grid of models and choosing the tolerances with the lowest average run time from those that conserved elements well and rarely failed. Despite this, there are no one-size-fits-all perfect tolerances and you may run into issues with different networks or models.
                        
        """)
        
        element_list = st.text_input('elements', value='H, N, C, O, S')
        element_list_arr = [_.strip() for _ in element_list.split(',')]
        conservation=uclchem.analysis.check_element_conservation(result_df, element_list=element_list_arr)
        st.write("Percentage change in total abundances:")
        st.dataframe(conservation)
    
    # with st.expander("Plotting result", expanded=True):
    with tab1:
        
        st.download_button(
            label="Download data as CSV",
            data=result_df.to_csv(),
            file_name='outputFile.csv',
            mime='text/csv',
        )
        
        st.markdown("""
            Note the use of $ symbols in the species list below, this gets the total ice abundance of a species. For two phase models, this is just the surface abudance but for three phase it is the sum of surface and bulk.
        """)
        
        species = st.text_input('species', value='H, H2, $H, $H2, H2O, $H2O, CO, $CO, $CH3OH, CH3OH')
        species_list = [_.strip() for _ in species.split(',')]
        
        # col1, col2 = st.columns(2)
        
        # xlim = col1.text_input('xlim', value='1e3, 1e6')
        # ylim = col2.text_input('ylim', value='1e-15, 1')
        
        abundances_dict = {}
        for specName in species_list:
            if specName[0] == "$":
                abundances_dict[specName] = result_df[specName.replace("$", "#")]
                if specName.replace("$", "@") in result_df.columns:
                    abundances_dict[specName] = abundances_dict[specName] + result_df[specName.replace("$", "@")]
            else:
                abundances_dict[specName] = result_df[specName]

        df = pd.DataFrame({'time': result_df["Time"]} | abundances_dict)
        df = df.set_index('time')
        fig = df.plot(
            title="Time (Year) vs Abundances", 
            log_y=True,
            labels=dict(index="Time / years", value="X<sub>{Species}</sub>", variable="Species"),
            # range_x=[float(_) for _ in xlim.split(',')], range_y=[float(_) for _ in ylim.split(',')]
        )
        st.plotly_chart(fig, use_container_width=True)

param_dict = {}
outSpecies: str = None
# file_saved = False
# status = -1

def main():
    
    global param_dict
    
    st.header("UCLCHEM")
    st.write("A Gas-grain Chemical Code for Clouds, Cores, and C-Shocks")
    st.divider()
    
    st.subheader("Fine-tuning the model")
    
    with st.expander("Parameters", expanded=True):
        parameters = get_parameters()
    with st.expander("Behavioural parameters"):
        behaviour_parameters = get_behavioural_parameters()
    # with st.expander("Input and Output parameters"):
    #     input_output_parameters = get_input_output_parameters()
        
    #     input_output_parameters_filtered = {
    #         key: set_loc(value)
    #         for key, value in input_output_parameters.items() if key != 'writeStep' and value
    #     }
    #     input_output_parameters_filtered['writeStep'] = input_output_parameters['writeStep']
        
        
    with st.expander("Integration Controls"):
        integration_controls = get_integration_controls()
    
    st.divider()
    
    st.markdown("""
        ## A Simple Cloud

        UCLCHEM's `cloud` model is a spherical cloud of isothermal gas. We can keep a constant density or have it increase over time following a freefall equation. This model is generally useful whenever you want to model a homogeneous cloud of gas under constant conditions. For example, in the inner parts of a molecular cloud where Av $\gtrsim$ 10 there are very few depth dependent processes. You may wish to model the whole of this UV shielded portion of the cloud with a single `cloud` model.

       _You can find a complete list of modifiable parameters and their default values as above._
        
    """)
    
    # species = st.text_input("outSpecies", value='SO, CO+')
    # out_species = [_.strip() for _ in species.split(',')]
    input_output_parameters = {
        "writeStep": 1,
        # "abundLoadFile": set_loc("abundance"),
        "abundSaveFile": set_loc("abundSaveFile"),
        # "columnFile": set_loc("columnFile"),
        "outputFile": set_loc("outputFile"),
    }
    param_dict = parameters | input_output_parameters | behaviour_parameters | integration_controls

    if not ('outputFile' in param_dict and param_dict['outputFile']):
        param_dict['outputFile'] = set_loc('outputFile')
    
    
    if st.button('Run calculations'):
        time_start = perf_counter()
        
        status, *_ = uclchem.model.cloud(param_dict=param_dict)
        if status < 0:
            return st.error('Error occured during calculation')
        st.success(f'Finished in {(perf_counter() - time_start):.2f} seconds')
        
        
    # st.button('Show result', on_click=simple_cloud_model_calc)
    simple_cloud_model_calc()
        
        
if __name__ == "__main__":
    
    main()
    about_page()
    