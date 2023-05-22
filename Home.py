import streamlit as st

# st.set_page_config(page_title="Home", page_icon=":rocket:", layout="wide")

"# A simple web interface to chemistry and astrochemistry tools"
"### Created and mainted by A.N. Marimuthu."
st.info("Choose a tool from the sidebar to get started")
# st.warning("This is a work in progress. Please report any bugs or issues to the author.")

with st.sidebar:
    """
        ## About this website
        
        I am Aravindh, a PhD student at the Radboud University, Nijmege, the Netherlands.
        I am a experimental spectroscopist of cold molecular ions that are relevant to interstellar medium (ISM).
        I am interested in astrochemistry and machine learning.
        
        This website provides a simple interface to some of the open-source tools and libraries in astrochemistry or chemistry in general, that are available in GitHub and other sources. 
        
        Check out the About section in sidebar on each page(s) for more information on the selected project.
        
        Check out my [Github profile](https://github.com/aravindhnivas)

    """
    
    st.warning("This is a work in progress. Please report any bugs or issues to the authour - aravindhnivas28+astrochemistry@gmail.com.")
    