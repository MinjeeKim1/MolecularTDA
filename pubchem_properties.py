import requests
import pandas as pd
from io import StringIO
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import time

def get_cid_from_name(drug_name):
    url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/cids/JSON"
    response = requests.get(url)

    if response.status_code == 200:
        data = response.json()
        cids = data.get("IdentifierList", {}).get("CID", [])
        if cids:
            return cids[0] 
    return None


def get_all_properties_to_excel(drug_names, output_filename="compounds_properties.xlsx"):
    properties = [
        "Title", "MolecularWeight", "XLogP", "ExactMass", "MonoisotopicMass",
        "TPSA", "HBondDonorCount", "HBondAcceptorCount", "RotatableBondCount",
        "HeavyAtomCount", "AtomStereoCount", "DefinedAtomStereoCount",
        "UndefinedAtomStereoCount", "BondStereoCount", "DefinedBondStereoCount",
        "UndefinedBondStereoCount", "CovalentUnitCount", "Volume3D",
        "XStericQuadrupole3D", "YStericQuadrupole3D", "ZStericQuadrupole3D",
        "FeatureAcceptorCount3D", "FeatureDonorCount3D", "FeatureAnionCount3D",
        "FeatureCationCount3D", "FeatureRingCount3D", "FeatureHydrophobeCount3D",
        "ConformerModelRMSD3D", "EffectiveRotorCount3D", "Fingerprint2D"
    ]

    all_data = []

    for drug in drug_names:
        print(f"Processing drug: {drug}")
        cid = get_cid_from_name(drug)

        if not cid:
            print(f"Could not find CID for drug: {drug}")
            continue

        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{','.join(properties)}/CSV"
        response = requests.get(url)

        if response.status_code == 200:
            csv_data = response.text
            data_frame = pd.read_csv(StringIO(csv_data))
            data_frame = data_frame.replace(r'^\s*$', 'N/A', regex=True)

            if "Disease" not in data_frame.columns:
                data_frame.insert(0, "Disease", "Huntington")

            all_data.append(data_frame.iloc[0])
        else:
            print(f" Error fetching properties for {drug} (CID {cid})")

    if all_data:
        combined_data_frame = pd.DataFrame(all_data)
        combined_data_frame.to_excel(output_filename, index=False)
        print(f" Saved all data to {output_filename}")
    else:
        print(" No data to save.")


def main(drug_names, output_filename="compounds_properties.xlsx"):
    get_all_properties_to_excel(drug_names, output_filename)

drug_names = ["Raxatrigine hydrochloride"]
if __name__ == "__main__":
    main(drug_names, output_filename="compounds_properties.xlsx")

