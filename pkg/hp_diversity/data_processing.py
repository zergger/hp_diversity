# diversity_analysis/data_processing.py
import pandas as pd
import numpy as np
from .calculations import d_chao, Cq, Uq, Sq, Vq

def load_data(file_path, sep="\t"):
    return pd.read_csv(file_path, sep=sep)

def calculate_alpha_diversity(df, d):
    Alpha = []
    grouped_data = df.groupby('ID')
    for group_name, group_data in grouped_data:
        group_array = np.nan_to_num(group_data.iloc[:, 1:].astype(float).values)
        result = d_chao(group_array, "alpha", d)
        Alpha.append([group_name] + [result])
    Alpha = np.array(Alpha)
    colnames = ["ID", f"q={d}"]
    Alpha = pd.DataFrame(Alpha, columns=colnames)
    return Alpha

def calculate_gamma_diversity(df, d):
    Gamma = []
    grouped_data = df.groupby('ID')
    for group_name, group_data in grouped_data:
        group_array = np.nan_to_num(group_data.iloc[:, 1:].astype(float).values)
        result = d_chao(group_array, "gamma", d)
        Gamma.append([group_name] + [result])
    Gamma = np.array(Gamma)
    colnames = ["ID", f"q={d}"]
    Gamma = pd.DataFrame(Gamma, columns=colnames)
    return Gamma

def calculate_beta_diversity(df1, df2, d):
    Res = pd.DataFrame(columns=["ID1", "ID2", "beta.D", "cq", "uq", "sq", "vq"])
    visited_pairs = set()
    grouped_data1 = df1.groupby('ID')
    grouped_data2 = df2.groupby('ID')
    for group_name1, group_data1 in grouped_data1:
        for group_name2, group_data2 in grouped_data2:
            if (group_name1, group_name2) in visited_pairs or (group_name2, group_name1) in visited_pairs:
                continue
            visited_pairs.add((group_name1, group_name2))            
            group_array1 = np.nan_to_num(group_data1.iloc[:, 1:].astype(float).values)
            group_array2 = np.nan_to_num(group_data2.iloc[:, 1:].astype(float).values)
            sda = np.vstack((group_array1, group_array2))
            beta = d_chao(sda, "beta", d)
            cqn = Cq(sda, beta, d, sda.shape[0])
            uqn = Uq(sda, beta, d, sda.shape[0]) 
            sqn = Sq(sda, beta, d, sda.shape[0]) 
            vqn = Vq(sda, beta, d, sda.shape[0])
            res = pd.DataFrame([[group_name1, group_name2, beta, cqn, uqn, sqn, vqn]], columns=["ID1", "ID2", "beta.D", "cq", "uq", "sq", "vq"])
            Res = pd.concat([Res, res], ignore_index=True)
    return Res
