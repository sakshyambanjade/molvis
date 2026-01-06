import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen
import streamlit.components.v1 as components


# Page Configuration
st.set_page_config(
    page_title="Molecule Viewer",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("Molecule Viewer")
st.markdown("Visualize molecular structures and basic properties from SMILES")


# Example molecules
EXAMPLE_MOLECULES = {
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Benzene": "c1ccccc1",
    "Ethanol": "CCO",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(=O)O1)O)O)O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}


# Sidebar controls
with st.sidebar:
    input_mode = st.radio(
        "Input Mode:",
        ["SMILES Text", "Draw Structure (Ketcher)"]
    )
    
    if input_mode == "SMILES Text":
        smiles_input = st.text_area(
            "SMILES:",
            placeholder="e.g., c1ccccc1",
            height=100
        )
    else:
        st.info("Draw your molecule below and click 'Get SMILES' to analyze")
        smiles_input = None


# JSME molecular editor integration
if input_mode == "Draw Structure (Ketcher)":
    st.subheader("Molecule Editor")
    
    jsme_html = """
    <!DOCTYPE html>
    <html>
    <head>
        <script src="https://peter-ertl.com/jsme/JSME_2022-05-01/jsme/jsme.nocache.js"></script>
    </head>
    <body>
        <div id="jsme_container"></div>
        <button onclick="getSmiles()" style="margin: 10px; padding: 10px 20px; background: #ff4b4b; color: white; border: none; border-radius: 5px; cursor: pointer; font-size: 14px;">
            Get SMILES from Drawing
        </button>
        <div id="smiles_output" style="margin: 10px; padding: 10px; background: #f0f0f0; border-radius: 5px; font-family: monospace;"></div>
        
        <script>
            var jsmeApplet;
            
            function jsmeOnLoad() {
                jsmeApplet = new JSApplet.JSME("jsme_container", "600px", "400px", {
                    "options": "query,hydrogens"
                });
            }
            
            function getSmiles() {
                if (jsmeApplet) {
                    var smiles = jsmeApplet.smiles();
                    document.getElementById('smiles_output').innerHTML = '<strong>SMILES:</strong> ' + smiles;
                    
                    // Copy to clipboard
                    navigator.clipboard.writeText(smiles).then(function() {
                        document.getElementById('smiles_output').innerHTML += ' <span style="color: green;">✓ Copied to clipboard!</span>';
                    });
                }
            }
        </script>
    </body>
    </html>
    """
    
    components.html(jsme_html, height=550)
    
    st.info("💡 Draw your molecule above, click 'Get SMILES from Drawing', then paste the SMILES into the text input on the left to analyze it.")

# Main analysis
if smiles_input and smiles_input.strip():
    # Parse and validate SMILES
    molecule = Chem.MolFromSmiles(smiles_input.strip())
    
    if molecule is not None:
        # Success message
        st.success(f"Valid SMILES: `{smiles_input.strip()}`")
        st.divider()
        
        # Display structure and properties
        col_structure, col_properties = st.columns(2)
        
        # Left Column - 2D Structure Visualization
        with col_structure:
            st.subheader("Molecular Structure")
            AllChem.Compute2DCoords(molecule)
            mol_image = Draw.MolToImage(molecule, size=(400, 400))
            st.image(mol_image, use_column_width=True)
        
        # Right Column - Key Properties
        with col_properties:
            st.subheader("Essential Properties")
            
            # Calculate properties once
            mol_weight = Descriptors.ExactMolWt(molecule)
            logp = Crippen.MolLogP(molecule)
            hbd = Descriptors.NumHDonors(molecule)
            hba = Descriptors.NumHAcceptors(molecule)
            tpsa = Descriptors.TPSA(molecule)
            rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
            aromatic_rings = Descriptors.NumAromaticRings(molecule)
            num_atoms = molecule.GetNumAtoms()
            
            # Display as metrics
            metric_col1, metric_col2 = st.columns(2)
            
            with metric_col1:
                st.metric("Molecular Weight (g/mol)", f"{mol_weight:.2f}")
                st.metric("H-Bond Donors", hbd)
                st.metric("Rotatable Bonds", rotatable_bonds)
            
            with metric_col2:
                st.metric("LogP (Lipophilicity)", f"{logp:.2f}")
                st.metric("H-Bond Acceptors", hba)
                st.metric("Aromatic Rings", aromatic_rings)
        
        st.divider()
        
        # Additional molecular descriptors
        st.subheader("Molecular Descriptors")
        
        desc_col1, desc_col2, desc_col3 = st.columns(3)
        
        with desc_col1:
            st.metric("Molar Refractivity", f"{Crippen.MolMR(molecule):.2f}")
            st.metric("Heavy Atoms", Descriptors.HeavyAtomCount(molecule))
        
        with desc_col2:
            st.metric("Ring Count", Descriptors.RingCount(molecule))
            st.metric("Heteroatoms", Descriptors.NumHeteroatoms(molecule))
        
        with desc_col3:
            st.metric("Aliphatic Rings", Descriptors.NumAliphaticRings(molecule))
            st.metric("Saturated Rings", Descriptors.NumSaturatedRings(molecule))
        
        st.divider()
        
        # 3D Structure Viewer
        with st.expander("3D Structure"):
            try:
                AllChem.EmbedMolecule(molecule, randomSeed=42)
                AllChem.MMFFOptimizeMolecule(molecule)
                
                # Generate SDF for 3D visualization
                sdf_block = Chem.MolToMolBlock(molecule)
                
                # Embed 3D viewer using py3Dmol
                html_code = f"""
                <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                <div id="viewer" style="width: 100%; height: 400px; position: relative;"></div>
                <script>
                    let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: 'white'}});
                    let pdbData = `{sdf_block}`;
                    viewer.addModel(pdbData, "sdf");
                    viewer.setStyle({{}}, {{cartoon: {{color: 'spectrum'}}}});
                    viewer.setStyle({{}}, {{stick: {{}}}});
                    viewer.zoomTo();
                    viewer.render();
                </script>
                """
                components.html(html_code, height=450)
            except:
                st.info("3D structure generation not available for this molecule")
    
    else:
        # Invalid SMILES
        st.error("Invalid SMILES string. Please check the syntax and try again.")
        st.markdown("""
        **Common SMILES issues:**
        - Missing ring closure numbers (e.g., `c1ccccc1` for benzene)
        - Unclosed parentheses
        - Invalid atom symbols
        """)

else:
    st.info("Enter a SMILES string in the sidebar to get started!")
    
# Footer
st.divider()