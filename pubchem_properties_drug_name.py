import requests
import pandas as pd
from io import StringIO
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import time

def get_cid_from_name(drug_name):
    try:
        import urllib.parse
        encoded_name = urllib.parse.quote(drug_name)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{encoded_name}/cids/JSON"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            data = response.json()
            cids = data.get("IdentifierList", {}).get("CID", [])
            if cids:
                return cids[0] 
        elif response.status_code == 404:
            print(f"  No compound found for: {drug_name}")
        else:
            print(f"  HTTP {response.status_code} for: {drug_name}")
    except requests.exceptions.RequestException as e:
        print(f"  Network error for {drug_name}: {e}")
    except Exception as e:
        print(f"  Unexpected error for {drug_name}: {e}")
    return None


def get_all_properties_to_excel(drug_names, output_filename="compounds_properties.xlsx"):
    properties = [
        "Fingerprint2D"
    ]

    all_data = []

    for drug in drug_names:
        print(f"Processing drug: {drug}")
        cid = get_cid_from_name(drug)

        if not cid:
            print(f"Could not find CID for drug: {drug}")
            continue

        try:
            url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/cid/{cid}/property/{','.join(properties)}/CSV"
            response = requests.get(url, timeout=10)

            if response.status_code == 200:
                csv_data = response.text
                data_frame = pd.read_csv(StringIO(csv_data))
                data_frame = data_frame.replace(r'^\s*$', 'N/A', regex=True)

                data_frame.insert(0, "DrugName", drug)

                if "Disease" not in data_frame.columns:
                    data_frame.insert(1, "Disease", "Huntington")

                all_data.append(data_frame.iloc[0])
                print(f"  Successfully processed {drug} (CID {cid})")
            else:
                print(f"  ✗ Error fetching properties for {drug} (CID {cid}): HTTP {response.status_code}")
        except requests.exceptions.RequestException as e:
            print(f"  Network error for {drug} (CID {cid}): {e}")
        except Exception as e:
            print(f"  Unexpected error for {drug} (CID {cid}): {e}")
        
        time.sleep(0.2)

    if all_data:
        combined_data_frame = pd.DataFrame(all_data)
        combined_data_frame.to_excel(output_filename, index=False)
        print(f"\nSaved {len(all_data)} compounds to {output_filename}")
        print(f"   Successfully processed: {len(all_data)}/{len(drug_names)} drugs")
    else:
        print("\nNo data to save - all compounds failed to process")


def main(drug_names, output_filename="compounds_properties.xlsx"):
    get_all_properties_to_excel(drug_names, output_filename)

df_pub = pd.read_csv("C:\\MolecularTDA\\getting_properties_scripts\\drugs.csv")
drug_names = df_pub['Drug Name']

if __name__ == "__main__":
    main(drug_names, output_filename="compounds_properties.xlsx")

