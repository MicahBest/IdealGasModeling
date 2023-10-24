import os
import pandas as pd

def txt_to_csv(file):
    """Converts txt file to csv"""
    data = []
    with open(file+".txt", 'r') as rf:
        for nrow, row in enumerate(rf):
            if nrow == 0:
                header = row.split()
            else:
                row_data = [float(s.replace(',', '')) for s in row.split()]
                row_dict = dict(zip(header, row_data))
                data.append(row_dict)

    df = pd.DataFrame(data)
    df.to_csv(file+".csv", index=False)

if __name__ == "__main__":
    path = "IGmodel/property-tables/"
    files = os.listdir(path)
    for file in files:
        if file.endswith(".txt"):
            file_path = os.path.join(path, file)
            txt_to_csv(file_path.replace(".txt", ''))

            # move txt to different directory
            new_path = path+"_txt/"
            txt_path = os.path.join(new_path, file)
            if not os.path.exists(new_path):
                os.mkdir(new_path)
            os.rename(file_path, txt_path)
