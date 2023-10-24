import os
import pandas as pd

def txt_to_csv(file, path="./", output_path="./"):
    """Converts txt file to csv"""
    data = []
    with open(path+file+".txt", 'r') as rf:
        for nrow, row in enumerate(rf):
            if nrow == 0:
                header = row.split()
            else:
                row_data = [float(s.replace(',', '')) for s in row.split()]
                row_dict = dict(zip(header, row_data))
                data.append(row_dict)

    df = pd.DataFrame(data)
    df.to_csv(output_path+file+".csv", index=False)

if __name__ == "__main__":
    input_path = "IGmodel/property-tables/_txt/"
    output_path = "IGmodel/property-tables/"
    files = os.listdir(input_path)
    for file in files:
        if file.endswith(".txt"):
            txt_to_csv(file.replace(".txt", ""), input_path, output_path)
