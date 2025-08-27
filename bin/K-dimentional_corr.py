#!/usr/bin/env python3
#  Copyright 2025 Jan Gorodkin, Sun Ying, Christian Anthon
# 
#  Use of this software is governed by the Business Source License
#  included in the file licenses/BSL.txt.
# 
#  As of the Change Date specified in that file, in accordance with
#  the Business Source License, use of this software will be governed
#  by the Apache License, Version 2.0, included in the file
#  licenses/APL.txt.

import sys
import numpy as np
import pandas as pd

def read_file(filename):
    df = pd.read_csv(filename, sep = '\t', header = 0)
    df = df[['editingeff', 'outcomefreq', 'predicted_editingeff', 'predicted_outcomefreq']]
    return df

def add_rank(df, columns):
    for column in columns:
        df.insert(df.shape[1], f'rank_{column}', df[column].rank())
        df.insert(df.shape[1], f'predicted_rank_{column}', df[f'predicted_{column}'].rank())
    return df

def calculate_covariance(x, y, k):
    xbar, ybar = x.mean(axis = 0), y.mean(axis = 0)
    return np.sum((x - xbar)*(y - ybar))/k

def calculate_rk_correlation(x, y, k):
    x, y = np.array(x), np.array(y)
    xx = calculate_covariance(x, x, k)
    yy = calculate_covariance(y, y, k)
    xy = calculate_covariance(x, y, k)
    corr = xy/np.sqrt(xx*yy)
    return corr

def get_R2_corr(df, columns, k):
    x = df[columns]
    y = df[[f'predicted_{column}' for column in columns]]
    corr = calculate_rk_correlation(x, y, k)
    return corr

def get_rho2_corr(df, columns, k):
    columns = [f'rank_{column}' for column in columns]
    x = df[columns]
    y = df[[f'predicted_{column}' for column in columns]]
    corr = calculate_rk_correlation(x, y, k)
    return corr

def main():
    filename = sys.argv[1]
    columns = ['editingeff', 'outcomefreq']
    k = len(columns)
    df = read_file(filename)
    df = add_rank(df, columns)
    R2 = get_R2_corr(df, columns, k)
    rho2 = get_rho2_corr(df, columns, k)
    print(f'R2:{R2:.4f}', f'rho2:{rho2:.4f}', sep = '\n')

if __name__ == '__main__':
    main()
