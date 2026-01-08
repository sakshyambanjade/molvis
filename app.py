import streamlit as st
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, Descriptors, Crippen
import streamlit.components.v1 as components


st.set_page_config(
    page_title="Molecule Viewer",
    layout="wide",
    initial_sidebar_state="expanded"
)

st.title("Molecule Viewer")
st.markdown("Visualize molecular structures and basic properties from SMILES")


EXAMPLE_MOLECULES = {
    "Caffeine": "CN1C=NC2=C1C(=O)N(C(=O)N2C)C",
    "Aspirin": "CC(=O)Oc1ccccc1C(=O)O",
    "Benzene": "c1ccccc1",
    "Ethanol": "CCO",
    "Glucose": "C([C@@H]1[C@H]([C@@H]([C@H](C(=O)O1)O)O)O)O",
    "Ibuprofen": "CC(C)Cc1ccc(cc1)C(C)C(=O)O",
}


with st.sidebar:
    st.header("Input")
    input_mode = st.radio(
        "Input Mode:",
        ["SMILES Text", "Draw Structure", "3D Viewer (PDB)"]
    )

    smiles_input = None
    pdb_id = None

    if input_mode == "SMILES Text":
        smiles_input = st.text_area(
            "SMILES:",
            placeholder="e.g., c1ccccc1",
            height=100,
        )

        if st.checkbox("Use example"):
            example_name = st.selectbox("Example molecule", list(EXAMPLE_MOLECULES.keys()))
            if example_name:
                smiles_input = EXAMPLE_MOLECULES[example_name]

    elif input_mode == "Draw Structure":
        st.info("Draw your molecule below and click 'Get SMILES' to analyze")

    elif input_mode == "3D Viewer (PDB)":
        pdb_id = st.text_input("PDB ID:", value="3PTB")


if input_mode == "Draw Structure":
    st.subheader("Kekule.js Molecule Editor")
    
    kekule_html = """
    <!DOCTYPE html>
    <html>
    <head>
        <link rel="stylesheet" type="text/css" href="https://unpkg.com/kekule/dist/themes/default/kekule.css" />
        <script src="https://unpkg.com/kekule/dist/kekule.min.js"></script>
        <style>
            body { 
                margin: 0; 
                padding: 20px; 
                font-family: sans-serif; 
            }
            #composer { 
                width: 100%; 
                height: 450px; 
                border: 2px solid #ddd; 
                border-radius: 8px; 
            }
            .button {
                margin-top: 15px;
                padding: 12px 24px;
                background: #ff4b4b;
                color: white;
                border: none;
                border-radius: 5px;
                cursor: pointer;
                font-size: 16px;
                font-weight: bold;
            }
            .button:hover {
                background: #e04242;
            }
            #output {
                margin-top: 15px;
                padding: 15px;
                background: #f0f2f6;
                border-radius: 5px;
                font-family: monospace;
                min-height: 30px;
                word-break: break-all;
            }
        </style>
    </head>
    <body>
        <div id="composer"></div>
        <button class="button" onclick="exportSmiles()">Get SMILES</button>
        <button class="button" onclick="clearEditor()" style="background: #666;">Clear</button>
        <div id="output"></div>
        
        <script>
            var composer;
            
            window.addEventListener('load', function() {
                composer = new Kekule.Editor.Composer(document.getElementById('composer'));
                composer.setDimension('100%', '450px');
                composer.setEnableOperHistory(true);
            });
            
            function exportSmiles() {
                var chemObj = composer.getChemObj();
                if (chemObj) {
                    try {
                        var smiles = Kekule.IO.saveFormatData(chemObj, 'smi');
                        smiles = smiles.trim();
                        
                        if (smiles && smiles !== '') {
                            document.getElementById('output').innerHTML = 
                                '<strong>SMILES:</strong> ' + smiles;
                            
                            navigator.clipboard.writeText(smiles).then(function() {
                                document.getElementById('output').innerHTML += 
                                    ' <span style="color: green;">✓ Copied to clipboard!</span>';
                            }).catch(function() {
                                document.getElementById('output').innerHTML += 
                                    ' <span style="color: orange;">(Copy manually)</span>';
                            });
                        } else {
                            document.getElementById('output').innerHTML = 
                                '<span style="color: red;">Please draw a molecule first!</span>';
                        }
                    } catch (e) {
                        document.getElementById('output').innerHTML = 
                            '<span style="color: red;">Error generating SMILES: ' + e.message + '</span>';
                    }
                } else {
                    document.getElementById('output').innerHTML = 
                        '<span style="color: red;">No molecule drawn!</span>';
                }
            }
            
            function clearEditor() {
                composer.setChemObj(null);
                document.getElementById('output').innerHTML = '';
            }
        </script>
    </body>
    </html>
    """
    
    components.html(kekule_html, height=650)
    
    st.info("Draw your molecule above, click 'Get SMILES', copy it, then switch to 'SMILES Text' mode and paste to analyze with RDKit.")


if input_mode == "3D Viewer (PDB)":
    st.subheader("3D Protein Viewer (Molstar)")

    if not pdb_id:
        pdb_id = "3PTB"

    molstar_html = f"""
    <!DOCTYPE html>
    <html>
    <head>
      <link rel="stylesheet" type="text/css" href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdbe-molstar-3.1.0.css">
      <style>
        body {{ margin: 0; padding: 0; }}
        #viewer-container {{ width: 100%; height: 500px; position: relative; }}
        #status {{ padding: 10px; text-align: center; color: #666; }}
      </style>
    </head>
    <body>
      <div id="status">Loading 3D structure...</div>
      <div id="viewer-container"></div>
      
      <script src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-plugin-3.1.0.js"></script>
      <script>
        const pdbId = "{pdb_id}".trim().toLowerCase() || "3ptb";
        
        const viewerInstance = new PDBeMolstarPlugin();
        
        const options = {{
          customData: {{
            url: `https://files.rcsb.org/download/${{pdbId}}.cif`,
            format: 'cif'
          }},
          alphafoldView: false,
          bgColor: {{ r: 255, g: 255, b: 255 }},
          hideControls: false,
          landscape: false,
          expanded: false,
          hideCanvasControls: ["selection", "animation"],
          sequencePanel: true
        }};
        
        viewerInstance.render(document.getElementById('viewer-container'), options)
          .then(() => {{
            document.getElementById('status').innerText = `Viewing: ${{pdbId.toUpperCase()}}`;
          }})
          .catch((error) => {{
            console.error('Molstar error:', error);
            document.getElementById('status').innerText = `Error loading ${{pdbId.toUpperCase()}}. Check if PDB ID is valid.`;
          }});
      </script>
    </body>
    </html>
    """

    components.html(molstar_html, height=600)
    st.info(f"Viewing PDB ID: {pdb_id.upper()} - Rotate with mouse, zoom with scroll. Try: 1CRN, 4HHB, 6LU7")


if input_mode == "SMILES Text":
    if smiles_input and smiles_input.strip():
        molecule = Chem.MolFromSmiles(smiles_input.strip())
        
        if molecule is not None:
            st.success(f"Valid SMILES: `{smiles_input.strip()}`")
            st.divider()

            col_structure, col_properties = st.columns(2)

            with col_structure:
                st.subheader("Molecular Structure")
                AllChem.Compute2DCoords(molecule)
                mol_image = Draw.MolToImage(molecule, size=(400, 400))
                st.image(mol_image, use_column_width=True)

            with col_properties:
                st.subheader("Essential Properties")

                mol_weight = Descriptors.ExactMolWt(molecule)
                logp = Crippen.MolLogP(molecule)
                hbd = Descriptors.NumHDonors(molecule)
                hba = Descriptors.NumHAcceptors(molecule)
                tpsa = Descriptors.TPSA(molecule)
                rotatable_bonds = Descriptors.NumRotatableBonds(molecule)
                aromatic_rings = Descriptors.NumAromaticRings(molecule)

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
            
            with st.expander("3D Structure"):
                try:
                    AllChem.EmbedMolecule(molecule, randomSeed=42)
                    AllChem.MMFFOptimizeMolecule(molecule)
                    sdf_block = Chem.MolToMolBlock(molecule)

                    html_code = f"""
                    <!DOCTYPE html>
                    <html>
                    <head>
                      <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
                    </head>
                    <body>
                      <div id="viewer" style="width: 100%; height: 400px; position: relative;"></div>
                      <script>
                        let viewer = $3Dmol.createViewer(document.getElementById("viewer"), {{backgroundColor: "white"}});
                        let sdfData = `{sdf_block}`;
                        viewer.addModel(sdfData, "sdf");
                        viewer.setStyle({{}}, {{stick: {{}}}});
                        viewer.zoomTo();
                        viewer.render();
                      </script>
                    </body>
                    </html>
                    """
                    components.html(html_code, height=450)
                except Exception:
                    st.info("3D structure generation not available for this molecule.")
        else:
            st.error("Invalid SMILES string. Please check the syntax and try again.")
            st.markdown("""
            Common SMILES issues:
            - Missing ring closure numbers (e.g., `c1ccccc1` for benzene)
            - Unclosed parentheses
            - Invalid atom symbols
            """)
    else:
        st.info("Enter a SMILES string in the sidebar to get started.")

st.divider()
