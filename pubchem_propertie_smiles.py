import requests
import pandas as pd
from io import StringIO
import time


def get_properties_from_smiles(smiles, disease="Huntington"):
    properties = ["Fingerprint2D"]

    try:
        import urllib.parse
        encoded_smiles = urllib.parse.quote(smiles)
        url = f"https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/smiles/{encoded_smiles}/property/{','.join(properties)}/CSV"
        response = requests.get(url, timeout=10)

        if response.status_code == 200:
            csv_data = response.text
            data_frame = pd.read_csv(StringIO(csv_data))
            data_frame = data_frame.replace(r'^\s*$', 'N/A', regex=True)

            data_frame.insert(0, "SMILES", smiles)

            if "Disease" not in data_frame.columns:
                data_frame.insert(1, "Disease", disease)

            return data_frame.iloc[0]
        elif response.status_code == 404:
            print(f"  No compound found for SMILES: {smiles}")
        else:
            print(f"  HTTP {response.status_code} for SMILES: {smiles}")
    except requests.exceptions.RequestException as e:
        print(f"  Network error for SMILES {smiles}: {e}")
    except Exception as e:
        print(f"  Unexpected error for SMILES {smiles}: {e}")

    return None


def get_all_properties_to_excel(smiles_list, output_filename="compounds_properties.xlsx"):
    all_data = []

    for smi in smiles_list:
        print(f"Processing SMILES: {smi}")
        row = get_properties_from_smiles(smi)

        if row is not None:
            all_data.append(row)
            print(f"  Successfully processed {smi}")
        else:
            print(f"  ✗ Failed to process {smi}")

        time.sleep(0.2)

    if all_data:
        combined_data_frame = pd.DataFrame(all_data)
        combined_data_frame.to_excel(output_filename, index=False)
        print(f"\nSaved {len(all_data)} compounds to {output_filename}")
        print(f"   Successfully processed: {len(all_data)}/{len(smiles_list)} SMILES")
    else:
        print("\nNo data to save - all SMILES failed to process")


def main(smiles_file, output_filename="compounds_properties.xlsx"):
    df_pub = pd.read_csv(smiles_file)
    smiles_list = df_pub['SMILES']
    get_all_properties_to_excel(smiles_list, output_filename)


if __name__ == "__main__":
    main("C:\\MolecularTDA\\getting_properties_scripts\\smiles.csv",
         output_filename="compounds_properties.xlsx")
