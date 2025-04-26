import csv
import math

def logp_normalize_and_average(csv_file_path):
    logp_values = []
    standard_value_index = None

    with open(csv_file_path, newline='') as csvfile:
        reader = csv.reader(csvfile, delimiter=';') 
        for row_num, row in enumerate(reader):
            if row_num == 0:
                try:
                    standard_value_index = row.index("Standard Value")
                    print(f"'Standard Value' column found at index {standard_value_index}")
                except ValueError:
                    print("Error: 'Standard Value' column not found in header.")
                    return
                continue

            if len(row) <= standard_value_index:
                continue

            raw_value = row[standard_value_index].strip()
            if not raw_value:
                continue

            try:
                value = float(raw_value)
                if value > 0:
                    logp = math.log10(value)
                    logp_values.append(logp)
            except ValueError:
                continue

    if logp_values:
        average_logp = sum(logp_values) / len(logp_values)
        print(f"\nAverage logP of 'Standard Value': {average_logp}")
    else:
        print("No valid 'Standard Value' entries found for logP computation.")

csv_file_path = "temp.csv"  
logp_normalize_and_average(csv_file_path)
