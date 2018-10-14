import pandas as pd

def read_dist(csv_file):
  df = pd.read_csv(csv_file)
  return df.values
  
def read_dist_dataframe(csv_file):
  df = pd.read_csv(csv_file)
  return df
  