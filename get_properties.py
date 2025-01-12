import requests
import pandas as pd
from io import StringIO
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.chrome.service import Service
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.options import Options
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
import time
from webdriver_manager.chrome import ChromeDriverManager
import xml.etree.ElementTree as ET
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen
import math
import re

from selenium import webdriver
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
from selenium.webdriver.common.by import By
from selenium.webdriver.support.ui import WebDriverWait
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.chrome.options import Options

def get_cid_from_name(drug_name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/Title/CSV'
    response = requests.get(url)
    
    if response.status_code == 200:
        try:
            data_frame = pd.read_csv(StringIO(response.text))
            if not data_frame.empty:
                cid = data_frame.iloc[0]["CID"]
                return cid
        except Exception as e:
            print(f"Error processing data for drug name {drug_name}: {e}")
            return None
    else:
        print(f"Error: Unable to fetch CID for {drug_name}. Status code: {response.status_code}")
        return None
    
def get_chembl_data(query):
    url = f'https://www.ebi.ac.uk/chembl/api/data/molecule/search?q={query}'
    print(url)
    response = requests.get(url)  
    try: 
        root = ET.fromstring(response.text) 
        molecule_properties = root.find(".//molecule/molecule_properties")
        apka = molecule_properties.find("cx_most_apka").text
        bpka = molecule_properties.find("cx_most_bpka").text
        return apka, bpka
    except Exception as e:
        print(f"Error processing data for {query}: {e}")
        return None, None
    

def get_name_from_cid(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/Title/CSV'
    response = requests.get(url)

    if response.status_code == 200:
        try:
            data_frame = pd.read_csv(StringIO(response.text))
            if not data_frame.empty:
                name = data_frame.iloc[0]["Title"]
                return name
        except Exception as e:
            print(f"Error processing data for CID {cid}: {e}")
            return None


def get_all_properties_to_excel(cids, output_filename="compounds_properties.xlsx"):
    properties = [
        "Title", #Name of the molecule
        "MolecularWeight", #Weight of the molecule
        "XLogP",  # Octanol-water partition coefficient
        "ExactMass",  # Precise molecular mass
        "MonoisotopicMass",  # Mass of the most abundant isotope
        "TPSA",  # Topological Polar Surface Area
        "HBondDonorCount",  # Hydrogen bond donors
        "HBondAcceptorCount",  # Hydrogen bond acceptors
        "RotatableBondCount",  # Bond flexibility
        "HeavyAtomCount",  # Molecular backbone size
        "AtomStereoCount",  # Stereochemistry (atoms)
        "DefinedAtomStereoCount",  # Defined stereochemistry (atoms)
        "UndefinedAtomStereoCount",  # Undefined stereochemistry (atoms)
        "BondStereoCount",  # Stereochemistry (bonds)
        "DefinedBondStereoCount",  # Defined stereochemistry (bonds)
        "UndefinedBondStereoCount",  # Undefined stereochemistry (bonds)
        "CovalentUnitCount",  # Number of covalent units
        "Volume3D",  # 3D spatial extent of the molecule
        "XStericQuadrupole3D",  # 3D shape/electron density (X component)
        "YStericQuadrupole3D",  # 3D shape/electron density (Y component)
        "ZStericQuadrupole3D",  # 3D shape/electron density (Z component)
        "FeatureAcceptorCount3D",  # 3D feature acceptors
        "FeatureDonorCount3D",  # 3D feature donors
        "FeatureAnionCount3D",  # 3D feature anions
        "FeatureCationCount3D",  # 3D feature cations
        "FeatureRingCount3D",  # 3D feature rings
        "FeatureHydrophobeCount3D",  # 3D feature hydrophobes
        "ConformerModelRMSD3D",  # Variation between conformers
        "EffectiveRotorCount3D",  # Molecular flexibility
        "Fingerprint2D"  # Encoded molecular substructures
    ]


    all_data = []

    for cid in cids:
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{','.join(properties)}/CSV"
        print(url)
        response = requests.get(url)

        if response.status_code == 200:
            csv_data = response.text
            data_frame = pd.read_csv(StringIO(csv_data))

            data_frame = data_frame.replace(r'^\s*$', 'N/A', regex=True)
            
            if "Disease" not in data_frame.columns:
                data_frame.insert(0, "Disease", "Huntington")  
            else:
                data_frame["Disease"] = "Huntington"  

            all_data.append(data_frame.iloc[0])  
            print(f"Data fetched for CID {cid}")
        else:
            print(f"Error: Unable to fetch data for CID {cid}. Status code: {response.status_code}")
        

        drugbank_data = fetch_drugbank_data(get_drugbank_id(cid))
        
        water_sol = drugbank_data.get("Water Solubility")
        if water_sol:
            all_data[-1]["WaterSolubility"] = water_sol
        else:
            all_data[-1]["WaterSolubility"] = "N/A"

        toxicity = drugbank_data.get("Rat acute toxicity")
        if toxicity:
            all_data[-1]["Toxicity"] = toxicity.split(" ")[0]
        else:
            all_data[-1]["Toxicity"] = "N/A"

        severity_grade = get_severity_grade(cid)
        if severity_grade:
            all_data[-1]["SeverityGrade"] = severity_grade
        else:
            all_data[-1]["SeverityGrade"] = "N/A"

        hepa_toxicity = get_hepa_toxicity(cid)
        if hepa_toxicity:
            all_data[-1]["Hepatotoxicity"] = hepa_toxicity
        else:
            all_data[-1]["Hepatotoxicity"] = "N/A"

        pka = get_chembl_data(get_name_from_cid(cid))
        pka_a = pka[0]
        if pka_a:
            all_data[-1]["APKA"] = pka[0]
        else:
            all_data[-1]["APKA"] = "N/A"
        pka_b = pka[1]
        if pka_b:
            all_data[-1]["BPKA"] = pka[1]
        else:
            all_data[-1]["BPKA"] = "N/A"


        aromatic_rings = find_aromatic_rings(cid)
        all_data[-1]["AromaticRings"] = aromatic_rings

        mol = Chem.MolFromSmiles(get_smiles(cid))
        atom_count = rdMolDescriptors.CalcNumAtoms(mol)
        all_data[-1]["AtomCount"] = atom_count

        caco_two = drugbank_data.get("Caco-2 permeable")
        if caco_two:
            all_data[-1]["Caco2Permeable"] = caco_two
        else:
            all_data[-1]["Caco2Permeable"] = "N/A"

        protein_binding_dictionary = drugbank_data.get("Protein binding")
        if protein_binding_dictionary:
            protein_binding_max = protein_binding_dictionary.get("max")
            protein_binding_min = protein_binding_dictionary.get("min")

            if protein_binding_max and protein_binding_min:
                all_data[-1]["ProteinBindingMin"] = protein_binding_min
                all_data[-1]["ProteinBindingMax"] = protein_binding_max
        else:
            all_data[-1]["ProteinBindingMin"] = "N/A"
            all_data[-1]["ProteinBindingMax"] = "N/A"
        
        blood_brain = drugbank_data.get("Blood Brain Barrier")
        if blood_brain:
            all_data[-1]["BloodBrainBarrier"] = blood_brain
        else:
            all_data[-1]["BloodBrainBarrier"] = "N/A"


        






    if all_data:
        combined_data_frame = pd.DataFrame(all_data)
        combined_data_frame.to_excel(output_filename, index=False)
        print(f"All data saved to {output_filename}")
    else:
        print("No data to save.")


def find_aromatic_rings(cid):
    smiles = get_smiles(cid)
    if not smiles: 
        print("Invalid SMILES string.")
        return []

    mol = Chem.MolFromSmiles(smiles)
    if not mol:
        print("Failed to convert SMILES string to molecule.")
        return []

    aromatic_rings = []
    ri = mol.GetRingInfo()
    for ring in ri.BondRings():
        if all(mol.GetBondWithIdx(bond_idx).GetIsAromatic() for bond_idx in ring):
            aromatic_rings.append([mol.GetBondWithIdx(b).GetBeginAtomIdx() for b in ring])
    return len(aromatic_rings)


def get_severity_grade(cid):
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--disable-gpu")
    options.add_argument("--start-maximized")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-extensions")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    driver.get(url)

    try:
        WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.XPATH, "//div[contains(text(), 'Severity Grade')]"))
        )
        severity_element = driver.find_element(
            By.XPATH, "//div[contains(text(), 'Severity Grade')]/following-sibling::div"
        )

        severity_grade = severity_element.text.strip()

        print(f"Severity Grade: {severity_grade}")
        return severity_grade
    except Exception as e:
        print(f"Error occurred while fetching Severity Grade: {e}")
        return None
    finally:
        driver.quit()

def get_hepa_toxicity(cid): 
    options = Options()
    options.add_argument("--headless")
    options.add_argument("--disable-gpu")
    options.add_argument("--start-maximized")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-extensions")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    driver.get(url)

    try:
        WebDriverWait(driver, 15).until(
            EC.presence_of_element_located((By.ID, "Hepatotoxicity"))
        )

        likelihood_element = driver.find_element(
            By.XPATH, "//section[@id='Hepatotoxicity']//p[contains(text(), 'Likelihood score')]"
        )

        likelihood_text = likelihood_element.text.strip()

        match = re.search(r"Likelihood score:\s*([A-Z])", likelihood_text)
        if match:
            likelihood_score = match.group(1)
            print(f"Likelihood Score (Letter): {likelihood_score}")
            return likelihood_score
        else:
            print("Likelihood score letter not found.")
            return None
    except Exception as e:
        print(f"Error occurred while fetching Likelihood Score: {e}")
        return None
    finally:
        driver.quit()


def get_smiles(cid):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/CanonicalSMILES/JSON'
    response = requests.get(url)
    
    if response.status_code == 200:
        data = response.json()
        if "PropertyTable" in data and "Properties" in data["PropertyTable"]:
            smiles_list = [item['CanonicalSMILES'] for item in data["PropertyTable"]["Properties"]]
            if smiles_list:
                return smiles_list[0] 
            else:
                return ""
        else:
            return ""
    else:
        return ""
    
def get_drugbank_id(cid):
    options = Options()
    options.add_argument("--headless") 
    options.add_argument("--disable-gpu")
    options.add_argument("--start-maximized")
    options.add_argument("--no-sandbox") 
    options.add_argument("--disable-extensions") 
    options.add_argument("--disable-plugins")  
    options.add_argument("--disable-images") 
    options.add_argument("--disable-javascript")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}"
    driver.get(url)

    try:
        drugbank_id_element = driver.find_element(By.XPATH, "//a[contains(@href, 'drugbank.ca')]")
        drugbank_id = drugbank_id_element.get_attribute("href")
        drugbank_id = drugbank_id.split("/")[-1]

        print(f"DrugBank ID: {drugbank_id}")
        return drugbank_id
    
    except Exception as e:
        print(f"Error occurred while fetching DrugBank ID: {e}")
    finally:
        driver.quit()
    
def calculate_esol(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return {"error": "Invalid SMILES"}

    mwt = Descriptors.MolWt(mol)

    clogp = Crippen.MolLogP(mol)

    rb = Descriptors.NumRotatableBonds(mol)
    
    num_aromatic_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetIsAromatic())
    num_heavy_atoms = Descriptors.HeavyAtomCount(mol)
    ap = num_aromatic_atoms / num_heavy_atoms if num_heavy_atoms > 0 else 0
    intercept = 0.16
    coeff_clogp = -0.63
    coeff_mwt = -0.0062
    coeff_rb = 0.066
    coeff_ap = -0.74

    logS = (intercept +
            coeff_clogp * clogp +
            coeff_mwt * mwt +
            coeff_rb * rb +
            coeff_ap * ap)
    
    solubility = 10 ** logS
    return solubility


def fetch_drugbank_data(drugbank_id):
    options = Options()
    options.add_argument("--disable-gpu")
    options.add_argument("--no-sandbox")
    options.add_argument("--disable-dev-shm-usage")
    options.add_argument("--disable-extensions")
    options.add_argument("--disable-plugins")
    options.add_argument("--disable-images")
    options.add_argument("--disable-javascript")
    options.add_argument("--headless")
    options.add_argument("user-agent=Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/114.0.0.0 Safari/537.36")

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)
    results = {}

    try:
        url = f"https://go.drugbank.com/drugs/{drugbank_id}"
        driver.get(url)

        WebDriverWait(driver, 20).until(
            EC.presence_of_element_located((By.ID, "drug-predicted-admet"))
        )

        # Fetching specific elements
        elements_to_fetch = {
            "Caco-2 permeable": "//td[text()='Caco-2 permeable']/following-sibling::td[2]",
            "Rat acute toxicity": "//td[text()='Rat acute toxicity']/following-sibling::td[1]",
            "hERG inhibition (predictor I)": "//td[text()='hERG inhibition (predictor I)']/following-sibling::td[2]",
            "hERG inhibition (predictor II)": "//td[text()='hERG inhibition (predictor II)']/following-sibling::td[2]",
            "Blood Brain Barrier": "//td[text()='Blood Brain Barrier']/following-sibling::td[2]",
        }

        for key, xpath in elements_to_fetch.items():
            try:
                element = WebDriverWait(driver, 20).until(
                    EC.presence_of_element_located((By.XPATH, xpath))
                )
                results[key] = element.text
            except Exception:
                results[key] = None  
                print(f"{key} not found on the page. Continuing to next fetch.")

        # Fetching protein binding data
        try:
            protein_binding_section = driver.find_element(
                By.XPATH, "//dt[text()='Protein binding']/following-sibling::dd"
            )

            protein_binding = protein_binding_section.text
            protein_binding = ''.join(filter(str.isdigit, protein_binding))
            min_max = {}

            if protein_binding:
                min_max["min"] = protein_binding
                min_max["max"] = protein_binding
            else:
                min_max["min"] = 'N/A'
                min_max["max"] = 'N/A'

            results["Protein binding"] = min_max
        except Exception as e:
            results["Protein binding"] = None
            print(f"Error occurred while fetching protein binding: {e}. Continuing.")

        # Fetching water solubility data
        try:
            water_solubility_xpath = '//td[text()="Water Solubility"]/following-sibling::td[1]'
            water_solubility_element = WebDriverWait(driver, 20).until(
                EC.presence_of_element_located((By.XPATH, water_solubility_xpath))
            )

            # Remove "mg/ml" and store only the value
            water_solubility = water_solubility_element.text.split(" ")[0]
            results['Water Solubility'] = water_solubility
        except Exception as e:
            results['Water Solubility'] = None
            print(f"Error occurred while fetching water solubility: {e}. Continuing.")

        # Handle case where water solubility is not found
        if results['Water Solubility'] is None:
            print("Water Solubility not found, dumping page source for debugging:")
            print(driver.page_source)

    except Exception as e:
        print(f"Error occurred while fetching data: {e}")
    finally:
        driver.quit()

    return results




def main(drug_names, output_filename="compounds_properties.xlsx"):
    cids = []
    for drug_name in drug_names:
        cid = get_cid_from_name(drug_name)
        if cid:
            cids.append(cid)
        else:
            print(f"Skipping {drug_name} due to missing CID.")
    
    if cids:
        get_all_properties_to_excel(cids, output_filename)
    else:
        print("No valid CIDs found. Exiting.")

#Replace with your drug names
drug_names = [
    "Amantadine", "Carbidopa", "Entacapone", "Levodopa", "Pramipexole", "Ramelteon", 
    "Ropinirole", "Tolcapone", "Alprazolam", "(3R,11bR)-3-(2-methylpropyl)-9,10-bis(trideuteriomethoxy)-1,3,4,6,7,11b-hexahydrobenzo[a]quinolizin-2-one", 
    "Dextromethorphan", "Digoxin", "Gsk-356278", "Haloperidol", "Idebenone", "(+)-Ketoconazole", 
    "Latrepirdine", "Nilotinib", "Omeprazole", "Risperidone", "Rolipram", "Sertraline", 
    "Ursodiol", "Warfarin", "Creatine", "Minocycline", "Acetylcysteine", "Atomoxetine", 
    "Bevantolol", "Biotin", "Bupropion", "Cannabidiol", "Citalopram", "Dronabinol", 
    "Fenofibrate", "Laquinimod", "Lithium Carbonate", "Mardepodect", "Mavoglurant", 
    "Memantine", "Neflamapimod", "N-[(E)-(4-bromo-3,5-dimethoxyphenyl)methylideneamino]-2-methoxy-2-(4-pyrazol-1-ylphenyl)acetamide", 
    "Phenylbutyric acid", "Pridopidine", "Srx246", "Thiamine", "Triheptanoin", "Ubidecarenone", 
    "Valproic Acid", "Quinidine", "Riluzole", "Tetrabenazine", "Tiapride", "Valbenazine"
]





main(drug_names)
