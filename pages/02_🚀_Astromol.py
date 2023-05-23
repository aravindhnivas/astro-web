import streamlit as st
import astromol as astro
import plotly.express as px
import pandas as pd

st.set_page_config(page_title="Astromol", layout="wide")

def intro():
    st.write(f"version {astro.__version__} (Last updated: {astro.updated()})")
    
    st.divider()
    
    """
    Molecules detected in space provide valuable insights into the physical and chemical processes that occur in the universe. These molecules can be detected using a variety of techniques, including radio telescopes and space probes.
    
    One of the most exciting aspects of detecting molecules in space is the potential for discovering new forms of life or habitable environments. For example, the discovery of water on Mars and the detection of organic molecules on comets and asteroids suggest that life may exist beyond Earth. Additionally, the detection of complex organic molecules in protoplanetary disks may provide clues as to how life evolved on Earth.

    Overall, the study of molecules in space is a rapidly growing field that has the potential to revolutionize our understanding of the universe and our place within it. As more advanced detection techniques are developed, we can expect to discover even more complex molecules and gain new insights into the origins and evolution of the cosmos.
    
    """
    
def telescope_summary(telescope: astro.Telescope):

    st.markdown(f"""
        - Short Name: {telescope.shortname}
        - Type: {telescope.type}
        - Diameter: {telescope.diameter}
        - Wavelengths: {', '.join(telescope.wavelength)}
        - Latitude: {telescope.latitude}
        - Longitude: {telescope.longitude}
        - Built: {telescope.built}
        - Decommissioned: {telescope.decommissioned}
    """)
    
    detects = [x.formula for x in astro.all_molecules if telescope in x.telescopes]
    st.write(f"{len(detects)} molecules detected With {telescope.name}")
    st.write(', '.join(detects))
    
def molecule_refs_summary(molecule: astro.Molecule):
    """
    Prints a nicely formatted list of references and notes to the terminal.
    """

    l_refs: list[str] = molecule.l_refs.split(";")
    d_refs: list[str] = molecule.d_refs.split(";")

    if molecule.notes != None:
        notes = molecule.notes.strip("*")

    st.markdown("**Detection Reference(s)**")
    strs = ""
    for x in range(len(d_refs)):
        strs += f"- {d_refs[x].strip()}\n"
    
    st.markdown(strs)

    st.markdown("**Laboratory Reference(s)**")
    
    strs = ""
    for x in range(len(l_refs)):
        strs += f"- {l_refs[x].strip()}\n"
    st.markdown(strs)
    
    if molecule.notes != None:
        st.write("_**Notes**_")
        st.write(notes)

    if molecule.isotopologues != None:
        iso_d_refs = molecule.isos_d_refs.split("[")
        del iso_d_refs[0]

        st.write("**Isotopologue Detection Reference(s)**")
        for x in iso_d_refs:
            st.write("- [" + x.strip())

    if molecule.ice == True or molecule.ice == "Tentative":
        st.write("**Ice Reference(s)**")
        # st.write("[Det] {}".format(molecule.ice_d_refs))
        # st.write("[Lab] {}".format(molecule.ice_l_refs))
        st.markdown(f"""
            - [Det] {molecule.ice_d_refs}
            - [Lab] {molecule.ice_l_refs}
        """)

    if molecule.ppd == True or molecule.ppd == "Tentative":
        st.write("**Protoplanetary Disks Reference(s)**")
        # st.write("[{}] {}".format(molecule.formula, molecule.ppd_d_refs))
        strs = f"- [{molecule.formula}] {molecule.ppd_d_refs}\n"
        
        if molecule.ppd_isos != None:
            for x in molecule.ppd_isos:
                # st.write(f"[{x.formula}] {x.ppd_d_refs}")
                strs += f"- [{x.formula}] {x.ppd_d_refs}\n"
        st.markdown(strs)

    if molecule.exgal == True or molecule.exgal == "Tentative":
        st.write("**External Galaxies Reference(s)**")
        st.write("- [{}] {}".format(molecule.formula, molecule.exgal_d_refs))

    if molecule.exo == True or molecule.exo == "Tentative":
        st.write("**Exoplanetary Atmospheres Reference(s)**")
        # st.write("[{}] {}".format(molecule.formula, molecule.exo_d_refs))
        strs = f"- [{molecule.formula}] {molecule.exo_d_refs}\n"
        if molecule.exo_isos != None:
            for x in molecule.exo_isos:
                # st.write(f"[{x.formula}] {x.exo_d_refs}")
                strs += f"- [{x.formula}] {x.exo_d_refs}\n"
        st.markdown(strs)
            
            
def molecule_summary(molecule):
    """
    Prints a summary of the information in the database for the molecule to the terminal.
    """
    
    attr_lists = ["Neutral", "Cation", "Anion", "Cyclic", "Radical"]
    attr_str = [attr_lists[ind] for ind, x in enumerate([molecule.neutral, molecule.cation, molecule.anion, molecule.cyclic, molecule.radical]) if x]
    
    st.markdown(f"""
        ### {molecule.name} ({molecule.formula})
        - Formula: {molecule.formula}
        - Name: {molecule.name}
        - Atoms: {molecule.natoms}
        - Mass: {molecule.mass} amu
        - Year Detected: {molecule.year}
        - Source(s): {', '.join([x.name for x in molecule.sources])}
        - Telescope(s) Used: {', '.join([x.shortname for x in molecule.telescopes])}
        - Attributes: {', '.join(attr_str)}
    """)
    
    other_envs = [molecule.ice, molecule.ppd, molecule.exgal, molecule.exo]
    
    other_str = ""
    if any(other_envs) == True:
        other_str = ""
        if molecule.ice == True:
            other_str += "Ices, "
        if molecule.ice == "Tentative":
            other_str += "Ices (Tentative), "
        if molecule.ppd == True:
            other_str += "Protoplanetary Disks, "
        if molecule.ppd == "Tentative":
            other_str += "Protoplanetary Disks (Tentative), "
        if molecule.exgal == True:
            other_str += "External Galaxies, "
        if molecule.exgal == "Tentative":
            other_str += "External Galaxies (Tentative), "
        if molecule.exo == True:
            other_str += "Exoplanetary Atmospheres, "
        if molecule.exo == "Tentative":
            other_str += "Exoplanetary Atmospheres (Tentative), "

        other_str = other_str.strip().strip(",")


    f"""
        {(f'- Known Isotopologues: {molecule.isotopologues}' if molecule.isotopologues else "")}
        {(f"- Also Detected In: {other_str}" if other_str else "")}
        {(f'- Source(s) of External Galaxy Detections: {molecule.exgal_sources}' if molecule.exgal_sources else "")}
        {(f'- Isotopologues Also Detected in Protoplanetary Disks: {",".join([x.formula for x in molecule.ppd_isos])}' if molecule.ppd_isos else "")}
        {(f'- Isotopologues Also Detected in Exoplanetary Atmospheres: {",".join([x.formula for x in molecule.exo_isos])}' if molecule.exo_isos else "")}
    """
    
    molecule_refs_summary(molecule)
    
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

def telescopes_full_lists():
    # full telescope table
    df = pd.DataFrame({
        'name': [_.name for _ in astro.all_telescopes],
        'shortname': [_.shortname for _ in astro.all_telescopes],
        'type': [_.type for _ in astro.all_telescopes],
        'diameter (in m)': [_.diameter for _ in astro.all_telescopes],
        'wavelength': [_.wavelength for _ in astro.all_telescopes],
        'latitude': [_.latitude for _ in astro.all_telescopes],
        'longitude': [_.longitude for _ in astro.all_telescopes],
        'built': [_.built for _ in astro.all_telescopes],
        'decommissioned': [_.decommissioned for _ in astro.all_telescopes],
        'ndetects': [_.ndetects for _ in astro.all_telescopes],
    })
    st.dataframe(df, use_container_width=True)
        
def main():
    intro()
    
    telescopes_tab, molecules_tab = st.tabs(["Telescopes", "Molecules"])
    
    total_molecules = len(astro.all_molecules)
    
    with telescopes_tab:
        
        st.markdown(
            f"""
            ### Total number of telescopes: {len(astro.all_telescopes)}
            
            - Detected {total_molecules} molecules in space.
            
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
        telescopes_full_lists()
        
        selected_telescope = st.selectbox("Select a telescope", [_.name for _ in astro.all_telescopes])
        selected_telescope_obj = [_ for _ in astro.all_telescopes if _.name == selected_telescope][0]
        telescope_summary(selected_telescope_obj)
        
    with molecules_tab:
        
        st.subheader(f"Total number of molecules detected in space: {total_molecules}")
        
        neutrals = [mol.name for mol in astro.all_molecules if mol.neutral]
        radicals = [mol.name for mol in astro.all_molecules if mol.radical]
        cations = [mol.name for mol in astro.all_molecules if mol.cation]
        cyclics = [mol.name for mol in astro.all_molecules if mol.cyclic]
        anions = [mol.name for mol in astro.all_molecules if mol.anion]
        fullerenes = [mol.name for mol in astro.all_molecules if mol.fullerene]
        pahs = [mol.name for mol in astro.all_molecules if mol.pah]
        
        df = pd.DataFrame({
            'neutrals': [len(neutrals), f"{(len(neutrals) / total_molecules * 100):.1f}"],
            'radicals': [len(radicals), f"{(len(radicals) / total_molecules * 100):.1f}"],
            'cations': [len(cations), f"{(len(cations) / total_molecules * 100):.1f}"],
            'cyclics': [len(cyclics), f"{(len(cyclics) / total_molecules * 100):.1f}"],
            'anions': [len(anions), f"{(len(anions) / total_molecules * 100):.1f}"],
            'fullerenes': [len(fullerenes), f"{(len(fullerenes) / total_molecules * 100):.1f}"],
            'PAHs': [len(pahs), f"{(len(pahs) / total_molecules * 100):.1f}"],
        }, index=['count', 'percentage (%)'], dtype=float).T
            
        st.dataframe(df, use_container_width=True)
        
        st.subheader("Molecule summary")
        selected_molecule = st.selectbox("Select a molecule", [_.name for _ in astro.all_molecules])
        selected_molecule_obj = [_ for _ in astro.all_molecules if _.name == selected_molecule][0]
        molecule_summary(selected_molecule_obj)
    
    about_page()

def about_page():
    st.divider()
    with st.sidebar:
        """
            ## About _astromol_
            

            A Database of Molecules Detected in Space
            
            The main python library is created and mainted by Dr. Brett A. McGuire
            
            [![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
            
            Code and models are available on [GitHub](https://github.com/bmcguir2/astromol), 
            and  [published paper](https://iopscience.iop.org/article/10.3847/1538-4365/aae5d2) 📃
            
            If you use `astromol` for your own work, please cite the Zenodo entry: [![DOI](https://zenodo.org/badge/360663606.svg)](https://zenodo.org/badge/latestdoi/360663606)
        """
        
if __name__ == "__main__":
    st.title("Astromol")
    # st.write("This web interface is created by A.N. Marimuthu.")
    main()
    