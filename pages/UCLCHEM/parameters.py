import streamlit as st
from pathlib import Path as pt

def get_parameters():
    
    # parameters
    st.markdown("### Parameters")
    
    col1, col2, col3, col4 = st.columns(4)
    
    initialTemp = col1.text_input("initialTemp", value="10", help="Initial gas temperature in Kelvin for all gas parcels in model")
    initialDens = col1.text_input("initialDens", value="1e2", help="Initial gas density in cm^-3 for all gas parcels in model")
    finalDens = col1.text_input("finalDens", value="1e5", help="Final gas density achieved via freefall")
    
    currentTime = col2.text_input("currentTime", value="0", help="Time at start of model in years")
    finalTime = col2.text_input("finalTime", value="5e6", help="Time to stop model in years, if not using `endAtFinalDensity` below")
    radfield = col2.text_input("radfield", value="1", help="Interstellar radiation field in Habing")
    
    zeta = col3.text_input("zeta", value="1", help="Cosmic ray ionisation rate as multiple of 1.3e-17 s^-1")
    rout = col3.text_input("rout", value="0.05", help="Outer radius of cloud being modelled in pc")
    rin = col3.text_input("rin", value="0", help="Minimum radial distance from cloud centre to consider")
    
    baseAv = col4.text_input("baseAv", value="2", help="Extinction at cloud edge, Av of a parcel at rout")
    points = col4.text_input("points", value="1", help="Number of gas parcels equally spaced between rin to rout to consider")
    
    return {
        "initialTemp": float(initialTemp),
        "initialDens": float(initialDens),
        "finalDens": float(finalDens),
        "currentTime": float(currentTime),
        "finalTime": float(finalTime),
        "radfield": float(radfield),
        "zeta": float(zeta),
        "rout": float(rout),
        "rin": float(rin),
        "baseAv": float(baseAv),
        "points": int(points)
    }
    

def get_behavioural_parameters():
    
    # parameters
    st.markdown("### Behavioural Controls")
    
    col1, col2, col3 = st.columns(3)
    
    h2desorb = col1.checkbox("h2desorb", value=True, help="Individually toggle non-thermal desorption due to H2 formation")
    crdesorb = col1.checkbox("crdesorb", value=True, help="Individually toggle non-thermal desorption due to cosmic rays")
    uvdesorb = col1.checkbox("uvdesorb", value=True, help="Individually toggle non-thermal desorption due to uv photons")
    desorb = col1.checkbox("desorb", value=True, help="Toggles all non-thermal desoprtion processes on or off")
    thermdesorb = col1.checkbox("thermdesorb", value=True, help="Toggle continuous thermal desorption")
    
    freefall = col2.checkbox("freefall", value=False, help="Controls whether models density increaes following freefall equation")
    endAtFinalDensity = col2.checkbox("endAtFinalDensity", value=False, help="Choose to end model at final density, otherwise end at final time")
    instantSublimation = col2.checkbox("instantSublimation", value=False, help="Toggle instantaneous sublimation of the ices at t")
    cosmicRayAttenuation = col2.checkbox("cosmicRayAttenuation", value=False, help="Use column density to attenuate cosmic ray ionisation rate following Padovani et al. 2018")
    improvedH2CRPDissociation = col2.checkbox("improvedH2CRPDissociation", value=False, help="Use H2 CRP dissociation rate from Padovani et al. 2018b")
    
    
    freezeFactor = col3.number_input("freezeFactor", value=1.0, help="Modify freefall rate by factor, usually to slow it")
    freefallFactor = col3.number_input("freefallFactor", value=1.0, help="Modify freeze out rate of gas parcels by this factor")
    ionModel = col3.text_input("ionModel", value="L", help="L/H model for cosmic ray attenuation Padovani et al. 2018")
    
    
    return {
        "freezeFactor": freezeFactor,
        "endAtFinalDensity": endAtFinalDensity,
        "freefall": freefall,
        "freefallFactor": freefallFactor,
        "desorb": desorb,
        "h2desorb": h2desorb,
        "crdesorb": crdesorb,
        "uvdesorb": uvdesorb,
        "thermdesorb": thermdesorb,
        "instantSublimation": instantSublimation,
        "cosmicRayAttenuation": cosmicRayAttenuation,
        "ionModel": ionModel,
        "improvedH2CRPDissociation": improvedH2CRPDissociation
    }
    
def get_input_output_parameters():
    # Input and Output parameters
    
    st.markdown("### Input and Output parameters")
    
    col1, col2 = st.columns(2)
    
    outputFile = col1.text_input("outputFile", value="outputFile", help="File to write full output of UCLCHEM. This includes physical parameter values and all abundances at every time step")
    columnFile = col1.text_input("columnFile", value="columnFile", help="File to write specific species abundances, see outSpecies")
    writeStep = col1.number_input("writeStep", value=1, help="Writing to columnFile only happens every writeStep timesteps")
    abundSaveFile = col2.text_input("abundSaveFile", value="abundSaveFile", help="File to store final abundances at the end of the model so future models can use them as the initial abundances. If not provided, no file will be produced")
    abundLoadFile = col2.text_input("abundLoadFile", value="", help="File from which to load initial abundances for the model, created through abundSaveFile. If not provided, the model starts from elemental gas")
    outSpecies = col2.text_input("outSpecies", value="", help="A space separated list of species to output to columnFile. Supplied as a separate list argument to most python functions")
    
    # loc = pt("./pages/UCLCHEM/outputs").absolute()
    # st.write(str(loc))
    # outputFile = str(loc / outputFile)
    # columnFile = str(loc / outputFile)
    
    # parameters = {"writeStep": writeStep}
    parameters = {
        "outputFile": outputFile,
        "columnFile": columnFile,
        "abundSaveFile": abundSaveFile,
        "abundLoadFile": abundLoadFile,
        "writeStep": int(writeStep),
        "outSpecies": outSpecies,
    }
        
    return parameters
    

def get_integration_controls():
    
    st.markdown("### Integration Controls")
    
    col1, col2, col3, col4 = st.columns(4)
    
    reltol = col1.text_input("reltol", value="1e-8", help="Relative tolerance for integration, see integration docs for advice")
    abstol_factor = col2.text_input("abstol_factor", value="1e-14", help="Absolute tolerance for integration is calculated by multiplying species abundance by this factor")
    abstol_min = col3.text_input("abstol_min", value="1e-25", help="Minimum value absolute tolerances can take")
    # MXSTEP = col4.number_input("MXSTEP", value=10000, help="Maximum steps allowed in integration before warning is thrown")
    
    return {
        "reltol": float(reltol),
        "abstol_factor": float(abstol_factor),
        "abstol_min": float(abstol_min),
        # "MXSTEP": int(MXSTEP)
    }
    
    
    