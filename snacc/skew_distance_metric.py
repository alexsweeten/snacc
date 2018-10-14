# A bijective function can be applied to the distance matrix generated using the default compression distance. For example, f(x) = -ln(1-x) transforms the usual range [0,1] to [0,+inf).

import numpy as np

from misc import read_dist_dataframe

def f_ln(x):
  return -np.log(1-x)
  
def f_inv(x):
  return x/(1-x)
  
def f_arctanh(x):
  return np.arctanh(x)

def main(csv_input, function, csv_output):
  D = read_dist_dataframe(csv_input)
  D_new = function(D)
  D_new.to_csv(csv_output, index=False)
  
def _test_skew():
  main("../test_dataset/distance_matrix_mysteryGenome1-8.csv", f_ln, "../test_dataset/distance_matrix_mysteryGenome1-8_skew.csv")
  
if __name__=="__main__":
  _test_skew()