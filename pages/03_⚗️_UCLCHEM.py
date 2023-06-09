import sys
import streamlit as st
import requests
from pages.UCLCHEM.parameters import (
    get_parameters, get_behavioural_parameters, 
    get_integration_controls
)
import pandas as pd

pd.options.plotting.backend = "plotly"
if len(sys.argv) > 1:
    devMode = sys.argv[1] == 'dev'
else:
    devMode = False
    
def about_page():
    
    st.divider()
    
    with st.sidebar:
        """
            ## About _UCLCHEM_ (v3.2.0)

            A Gas-Grain Chemical Code for astrochemical modelling
            
            UCLCHEM is a gas-grain chemical code for astrochemical modelling. It propagates the abundances of chemical species through a network of user-defined reactions according to the physical conditions of the gas.
            
            visit [here](https://uclchem.github.io/) for more information.
            
            [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
            
            Code and models are available on [GitHub](https://github.com/uclchem/UCLCHEM), 
            and  [published paper](https://iopscience.iop.org/article/10.3847/1538-3881/aa773f)
            📃
            
            Teams at UCL and Leiden Observatory are working with international collaborators to extend and improve UCLCHEM. Check our their [blog](https://uclchem.github.io/blog/index.html) for recent updates.
            
        """


def compute_data(api='api/simple_model'):
    print(f"Running API in {'development' if devMode else 'production'} mode")
    if devMode:
        URL = "http://localhost:9090"
    else:
        URL = st.secrets['UCLCHEM_API_URL']
        
    response = requests.post(f'{URL}/{api}', json=param_dict, headers={'Content-Type': 'application/json'})
    
    if response.status_code == 200:
        response_json = response.json()
        if response_json['status'] > 0:
            st.success(f'Calculation finished in {response_json["time"]:.2f} seconds')
            return response_json['results']
    else:
        response.raise_for_status()
        return None


def simple_cloud_model_calc():
    
    global param_dict
    
    st.markdown("""
        ## A Simple Cloud

        UCLCHEM's `cloud` model is a spherical cloud of isothermal gas. We can keep a constant density or have it increase over time following a freefall equation. This model is generally useful whenever you want to model a homogeneous cloud of gas under constant conditions. For example, in the inner parts of a molecular cloud where Av $\gtrsim$ 10 there are very few depth dependent processes. You may wish to model the whole of this UV shielded portion of the cloud with a single `cloud` model.

       _You can find a complete list of modifiable parameters and their default values as above._
        
    """)
    
    input_output_parameters = {
        "writeStep": 1
    }
    
    param_dict = param_dict | input_output_parameters
    if st.button('Run calculations'):
        st.session_state['simple_model_results'] = compute_data(api='api/simple_model')
    
    simple_model_results_page()
        

def simple_model_results_page():
            
    if 'simple_model_results' not in st.session_state:
        st.session_state['simple_model_results'] = None
        return
    
    result = st.session_state['simple_model_results']
    if not result: return st.warning('No results to display')
    tab1, tab2 = st.tabs(['Results', 'Elemental conservation'])
    
    with tab1:
        result_df = pd.read_json(result['full_output'])
        st.download_button(
            label="Download data as CSV",
            data=result_df.to_csv(),
            file_name='outputFile.csv',
            mime='text/csv',
        )
        
        st.markdown("""
            ### Species
            Note the use of $ symbols in the species list below, this gets the total ice abundance of a species. For two phase models, this is just the surface abudance but for three phase it is the sum of surface and bulk.
        """)
    
        species = st.text_input('Enter species', value='H, H2, $H, $H2, H2O, $H2O, CO, $CO, $CH3OH, CH3OH')
        species_list = [_.strip() for _ in species.split(',')]
          
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
        )
        st.plotly_chart(fig, use_container_width=True)
    
    with tab2:
        st.markdown("""
            To test whether we conserve elements. Each entry gives the change in the total abundance of an element as a percentage of the original abundance. In an ideal case, these values are 0\% indicating the total abundance at the end of the model is exactly the same as the total at the start.
            
            Changes of less than 1\% are fine for many cases but if they are too high, you could consider changing the `reltol` and `abstol` parameters that control the integrator accuracy. They are error tolerance so smaller values lead to smaller errors and (usually) longer integration times. The default values were chosen by running a large grid of models and choosing the tolerances with the lowest average run time from those that conserved elements well and rarely failed. Despite this, there are no one-size-fits-all perfect tolerances and you may run into issues with different networks or models.
                        
        """)
        
        conservation = result['conservation']
        st.write("Percentage change in total abundances:")
        st.dataframe(conservation)
    
    
param_dict = {}
outSpecies: str = None


def main():
    
    global param_dict
    
    st.header("UCLCHEM v3.2.0")
    st.write("A Gas-grain Chemical Code for Clouds, Cores, and C-Shocks")
    st.divider()
    st.subheader("Fine-tuning the model")
    
    with st.expander("Parameters", expanded=True):
        parameters = get_parameters()
    with st.expander("Behavioural parameters"):
        behaviour_parameters = get_behavioural_parameters()
    
    with st.expander("Integration Controls"):
        integration_controls = get_integration_controls()
        
    param_dict = parameters | behaviour_parameters | integration_controls
    st.divider()
    
    tab1, = st.tabs(['Simple cloud model'])
    with tab1:
        simple_cloud_model_calc()
    
if __name__ == "__main__":
    
    main()
    about_page()
    