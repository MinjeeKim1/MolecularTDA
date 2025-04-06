import csv
import requests
from rdkit import Chem
from selenium import webdriver
from selenium.webdriver.common.by import By
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.chrome.service import Service
from webdriver_manager.chrome import ChromeDriverManager
import time
import os
import pandas as pd

csv_file_path = 'drugs.csv'
drug_names = []

# get drug names from csv
with open(csv_file_path, mode='r', newline='', encoding='utf-8') as csvfile:
    reader = csv.DictReader(csvfile)
    for row in reader:
        if 'Drug Name' in row:
            drug_names.append(row['Drug Name'])

# get smiles from pubchem
def get_smiles_from_name(drug_name):
    url = f'https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{drug_name}/property/IsomericSMILES/JSON'
    response = requests.get(url)
    if response.status_code == 200:
        try:
            return response.json()['PropertyTable']['Properties'][0]['IsomericSMILES']
        except (KeyError, IndexError):
            return None
    else:
        return None

# get smiles for each drug name
smiles_list = []
for name in drug_names:
    smiles = get_smiles_from_name(name)
    smiles_list.append(smiles)

#convert smiles to mol
mol_list = []
for smiles in smiles_list:
    if smiles:
        mol = Chem.MolFromSmiles(smiles)
        mol_list.append(mol)
    else:
        mol_list.append(None)

print(smiles_list[0])

# get properties using http://www.swissadme.ch/index.php webscrapping
# Function to get properties from SwissADME
def get_properties_from_swissadme(smiles, name):
    download_dir = os.path.abspath("downloads")

    options = webdriver.ChromeOptions()
    prefs = {
        "download.default_directory": download_dir,
        "download.prompt_for_download": False,
        "directory_upgrade": True,
        "safebrowsing.enabled": True
    }
    options.add_experimental_option("prefs", prefs)

    driver = webdriver.Chrome(service=Service(ChromeDriverManager().install()), options=options)

    #open the web page
    driver.get("http://www.swissadme.ch/index.php")
    time.sleep(2)

    #find the smiles input box
    smiles_input = driver.find_element(By.ID, "smiles")
    smiles_input.clear()
    smiles_input.send_keys(f"{smiles}") 
    print("SMILES entered")

    #pres the run btn
    run_button = driver.find_element(By.ID, "submitButton")
    run_button.click()
    print("Run button clicked")
    time.sleep(20)

    # clicking the download csv button
    download_button = driver.find_element(By.XPATH, "//a[contains(@href, 'results/')]")
    download_button.click()
    time.sleep(5)

    # renaming the downloaded csv file
    downloaded_file = os.path.join(download_dir, "swissadme.csv")
    renamed_file = os.path.join(download_dir, f"{name}_properties.csv")
    if os.path.exists(downloaded_file):
        os.rename(downloaded_file, renamed_file)

    print(f"csv file download")

    driver.quit()

#run for 1 smiles
get_properties_from_swissadme(smiles_list[0], drug_names[0])
print("Properties downloaded")





