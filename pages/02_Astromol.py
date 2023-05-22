import streamlit as st
import astromol as astro
import plotly.express as px

st.set_page_config(page_title="Astromol", layout="wide")

def intro():
    st.write(f"version {astro.__version__} (Last updated: {astro.updated()})")
    
    st.divider()
    
    """
    Molecules detected in space provide valuable insights into the physical and chemical processes that occur in the universe. These molecules can be detected using a variety of techniques, including radio telescopes and space probes.
    
    One of the most exciting aspects of detecting molecules in space is the potential for discovering new forms of life or habitable environments. For example, the discovery of water on Mars and the detection of organic molecules on comets and asteroids suggest that life may exist beyond Earth. Additionally, the detection of complex organic molecules in protoplanetary disks may provide clues as to how life evolved on Earth.

    Overall, the study of molecules in space is a rapidly growing field that has the potential to revolutionize our understanding of the universe and our place within it. As more advanced detection techniques are developed, we can expect to discover even more complex molecules and gain new insights into the origins and evolution of the cosmos.
    
    """
    
def telescope_summary(telescope):

    # st.write('-'*len(telescope.name))
    # st.write(f"{telescope.name}")
    # st.write('-'*len(telescope.name))

    st.write(f"{'Type':11}\t{telescope.type}")
    if telescope.diameter:
        st.write(f"{'Diameter':11}\t{telescope.diameter}")
    st.write(f"{'Wavelengths':11}\t{', '.join(telescope.wavelength)}")
    if telescope.latitude:
        st.write(f"{'Latitude':11}\t{telescope.latitude}")
    if telescope.longitude:
        st.write(f"{'Longitude':11}\t{telescope.longitude}")
    st.write(f"{'Built':11}\t{telescope.built}")
    if telescope.decommissioned:
        st.write(f"{'Decommissioned':11}\t{telescope.decommissioned}")

    # mol_str = f"Molecules Detected With {telescope.name}"
    # st.write("\n" + "-"*len(mol_str))
    # st.write(mol_str)
    # st.write("-"*len(mol_str))
    

    detects = [x.formula for x in astro.all_molecules if telescope in x.telescopes]
    st.write(f"{len(detects)} molecules detected With {telescope.name}")
    st.write(', '.join(detects))
    
def get_telescopes():
    mainKey = f'updated: {astro.updated()}'
    data = {
        'names': [mainKey],
        'parents': [''],
        'values': [0],
    }

    data['names'] += ['Telescopes']
    data['parents'] += [mainKey]
    data['values'] += [1]

    for telescope in astro.all_telescopes:
        
        if not telescope.name:
            continue
            
        # add telescope
        data['names'].append(telescope.name)
        data['parents'].append('Telescopes')
        # data['parents'].append(telescope.wavelength)
        data['values'].append(telescope.ndetects)
        
        # add observed sources in that telescope
        observed_molecules = []
        for molecule in telescope.mols:
            if not molecule.name:
                continue
            
            molecule_name = f"{molecule.formula} ({molecule.name})"
            observed_molecules.append(molecule_name)
        
        # print(f"{telescope.name}: {len(observed_molecules)=} {len(telescope.mols)=}")
        if not observed_molecules:
            continue
        data['names'] += observed_molecules
        data['parents'] += len(observed_molecules)*[telescope.name]
        value = len(observed_molecules) / 360 * 100
        data['values'] += len(observed_molecules)*[value]

    fig = px.sunburst(
        data, **data,
        width = 1500,
        height = 1000
    )
    
    return fig


def main():
    intro()
    
    telescopes_tab, molecular_catagories_tab = st.tabs(["Telescopes", "Molecular catagories"])
    
    with telescopes_tab:
        # st.subheader("Telescopes")
        st.markdown(
            f"""
            ### Total number of telescopes: {len(astro.all_telescopes)}
            
            - Detected {len(astro.all_molecules)} molecules in space.
            
            The following shows the number of telescopes in each wavelength range.
            
            
            - IR : {len([_ for _ in astro.all_telescopes if "IR" in _.wavelength])}
            - _cm_ : {len([_ for _ in astro.all_telescopes if "cm" in _.wavelength])}
            - _mm_ : {len([_ for _ in astro.all_telescopes if "mm" in _.wavelength])}
            - _submm_ : {len([_ for _ in astro.all_telescopes if "sub-mm" in _.wavelength])}
            - _UV_ : {len([_ for _ in astro.all_telescopes if "UV" in _.wavelength])}
            - _Vis_ : {len([_ for _ in astro.all_telescopes if "Vis" in _.wavelength])}
            
            _NOTE_: Some telescope have multiple wavelength ranges.
            
            The following figure shows the number of molecules detected by each telescope. The size of the circle is proportional to the number of molecules detected by that telescope. The larger the circle, the more molecules have been detected by that telescope.
            """
        )
        st.plotly_chart(get_telescopes(), use_container_width=True)
        
        st.subheader("Telescopes summary")
        selected_telescope = st.selectbox("Select a telescope", [_.name for _ in astro.all_telescopes])
        selected_telescope_obj = [_ for _ in astro.all_telescopes if _.name == selected_telescope][0]
        # st.subheader(f"{selected_telescope_obj.name} summary")
        telescope_summary(selected_telescope_obj)
        
    with molecular_catagories_tab:
        st.subheader("Molecular catagories")
        
    
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
    