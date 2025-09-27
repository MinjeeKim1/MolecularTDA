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

    driver.get("http://www.swissadme.ch/index.php")
    time.sleep(2)

    smiles_input = driver.find_element(By.ID, "smiles")
    smiles_input.clear()
    smiles_input.send_keys(f"{smiles}") 
    print("SMILES entered")

    run_button = driver.find_element(By.ID, "submitButton")
    run_button.click()
    print("Run button clicked")
    time.sleep(20)

    download_button = driver.find_element(By.XPATH, "//a[contains(@href, 'results/')]")
    download_button.click()
    time.sleep(5)

    downloaded_file = os.path.join(download_dir, "swissadme.csv")
    renamed_file = os.path.join(download_dir, f"{name}_properties.csv")
    if os.path.exists(downloaded_file):
        os.rename(downloaded_file, renamed_file)

    print(f"csv file download")

    driver.quit()

def main(input_csv):
    df = pd.read_csv(input_csv)
    smiles_list = df['Obabel'].tolist()
    names_list = df['Drug Name'].tolist()

    for smiles, name in zip(smiles_list, names_list):
        print(f"Processing {name} with SMILES: {smiles}")
        get_properties_from_swissadme(smiles, name)
        time.sleep(5) 

if __name__ == "__main__":
    input_csv = "i.csv"  
    main(input_csv)
