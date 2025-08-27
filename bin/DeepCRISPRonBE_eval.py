#!/usr/bin/env python3
#  Copyright 2025 Jan Gorodkin, Sun Ying, Christian Anthon, Giulia Corsi
# 
#  Use of this software is governed by the Business Source License
#  included in the file licenses/BSL.txt.
# 
#  As of the Change Date specified in that file, in accordance with
#  the Business Source License, use of this software will be governed
#  by the Apache License, Version 2.0, included in the file
#  licenses/APL.txt.


import re
import sys
import glob
import numpy as np
import pandas as pd
from tensorflow import keras
from itertools import combinations

def check_weight_is_number(weight, weight_error):
    tmp = weight
    tmp = tmp.replace('-', '').replace('.', '')
    if not tmp.isdigit():
        print(weight_error, sep = '\n', file=sys.stderr)
        print(f'Current {weight} includes abnormal characters', sep = '\n', file=sys.stderr)
        sys.exit(1)

def check_the_number_of_weight(weight, editor, weight_error):
    weight = [float(i) for i in weight.split('-')]
    n = len(weight)
    if (editor == 'ABE' and n != 5) or (editor == 'CBE' and n != 3):
        print(weight_error, sep = '\n', file=sys.stderr)
        print(f'Current {weight} does not include correct number of weights', sep = '\n', file=sys.stderr)
        sys.exit(1)
    return weight

def check_sum_of_weight(weight, weight_error):
    s = sum(weight)
    if s<0.99 or s>1:
        print(weight_error, sep = '\n', file=sys.stderr)
        print(f'Current sum of {weight} is {s} not 1', sep = '\n', file=sys.stderr)
        sys.exit(1)

def check_weight(editor, weight):
    weight_error = '\n'.join(['Please provide the correct format for weight assignment.',
                       '5 numbers for ABE or 3 numbers for CBE are required.',
                       'These numbers are connected by "-" and the sum of the weights should be 1.',
                       'for example: "0.5-0.5-0-0-0" for ABE and "0.5-0.5-0" for CBE.'])
    check_weight_is_number(weight, weight_error)
    weight = check_the_number_of_weight(weight, editor, weight_error)
    check_sum_of_weight(weight, weight_error)
    return weight

def get_setting(editor):
    if editor not in ['ABE', 'CBE']:
        print(f'Currently {editor} prediction is not supported by our model', 
                'Please choose ABE or CBE for efficiency prediction', sep = '\n', file=sys.stderr)
        sys.exit(1)
    edited = {'ABE': 'G', 'CBE': 'T'}
    seq_length = 30
    editing_window = [4+3-1, 4+10]
    target_nucleotide = editor[0]
    edit_nucleotide = edited[editor]
    if len(sys.argv) > 7:
        weight = sys.argv[7]
        weight = check_weight(editor, weight)
    else:
        if editor == 'ABE':
            weight = [0.5, 0.5, 0, 0, 0]
        else:
            weight = [0.5, 0.5, 0]
    return seq_length, editing_window, target_nucleotide, edit_nucleotide, weight

def running_log(path_to_out):
    output = f'{path_to_out}/log.txt'
    w = open(output, 'w')
    return w

def read_fasta_file(filename):
    seqs, num = {}, 0
    with open(filename) as f:
        for line in f:
            line = line.rstrip()
            if line[0] == '>':
                seq_id = line[1:]
                num += 1
            else:
                seq = line
                seqs[seq_id] = seq
    return seqs, num

def check_unique(seqs, num):
    if len(seqs) < num:
        print('sequence identifiers must be unique', file=sys.stderr)
        sys.exit(1)

def check_length(seqs, seq_length, remove_targets, w):
    tmp = []
    print(f'#check all target sequences are 30 nt', file = w)
    for seq_id, seq in seqs.items():
        if len(seq) != seq_length:
            print(f'{seq_id} is not a sequence of length {seq_length}', file = w)
            tmp.append(seq_id)
    if len(tmp) > 0:
        remove_targets += tmp
    else:
        print(f'#pass', file = w)
    print(file = w)

def check_PAM(seqs, seq_length, remove_targets, w):
    tmp = []
    print('#check all sequences with PAM NGG', file = w)
    for seq_id, seq in seqs.items():
        if seq[seq_length-5: seq_length-3] != 'GG':
            print(f"{seq_id}'s PAM is not NGG", file = w)
            tmp.append(seq_id)
    if len(tmp) > 0:
        remove_targets += tmp
    else:
        print('#pass', file = w)
    print(file = w)

def check_edtiable_nucleotide(seqs, editor, editing_window, remove_targets, w):
    tmp = []
    target_nucleotide = editor[0]
    print(f'#check all sequences with at least one {target_nucleotide} in editing window from position 3 to 10 in gRNAs', file = w)
    for seq_id, seq in seqs.items():
        n = seq[editing_window[0]: editing_window[1]].count(target_nucleotide)
        if n == 0:
            print(f'{seq} without {target_nucleotide} in editing window(position 3-10 in gRNAs)', file = w)
            tmp.append(seq_id)
    if len(tmp) > 0:
        remove_targets += tmp
    else:
        print('#pass', file = w)
    print(file = w)

def remove_seqs(seqs, remove_targets):
    new_seqs = {seq_id: seq for seq_id, seq in seqs.items() if seq_id not in remove_targets}
    return new_seqs

def check_data_available(seqs):
    if len(seqs) == 0:
        print(
            '',
            'Error !!!  No sequences for prediction',
            'CRISPRon-BE has the following requirements for the input:',
            '30nt target DNA sequence (4nt+20nt protospacer+3nt PAM+3nt)',
            'PAM should be NGG motif',
            'At least one editable nucleotide (A for ABE, C for CBE) in the editing window (position 3 to 10 in PAM-distal end)',
            '',
            sep = '\n', file=sys.stderr)
        sys.exit(1)

def obtain_seqs(editor, filename, seq_length, editing_window, w):
    remove_targets = []
    seqs, num = read_fasta_file(filename)
    check_unique(seqs, num)
    check_length(seqs, seq_length, remove_targets, w)
    check_PAM(seqs, seq_length, remove_targets, w)
    check_edtiable_nucleotide(seqs, editor, editing_window, remove_targets, w)
    seqs = remove_seqs(seqs, remove_targets)
    check_data_available(seqs)
    return seqs

def get_seq_edited_positions(seq, target_nucleotide, editing_window):
    raw_locations = [i.start() for i in re.finditer(target_nucleotide, seq[editing_window[0]:editing_window[1]])]
    locations = [i+editing_window[0] for i in raw_locations]
    return raw_locations, locations

def get_specific_outcome(seq, array, edit_nucleotide):
    outcome = [edit_nucleotide if i in array else n for i,n in enumerate(list(seq))]
    outcome = ''.join(outcome)
    return outcome

def generate_one_target_outcome_feature(outcome, raw_locations, edit_nucleotide, editing_window):
    outcome_feature = [1 if i in raw_locations and n==edit_nucleotide else 0 for i,n in enumerate(list(outcome[editing_window[0]:editing_window[1]]))]
    return outcome_feature

def two_key_dic(dic, k1, k2, v):
    if not k1 in dic:
        dic[k1] = {}
    if not k2 in dic[k1]:
        dic[k1][k2] = v

def generate_one_target_outcomes(seq, outcome_features, target_nucleotide, edit_nucleotide, editing_window):
    outcomes = []
    raw_locations, locations = get_seq_edited_positions(seq, target_nucleotide, editing_window)
    location_combs = sum([list(map(list, combinations(locations, i))) for i in range(len(locations) + 1)], [])
    for array in location_combs:
        outcome = get_specific_outcome(seq, array, edit_nucleotide)
        if outcome == seq:
            continue
        outcome_feature = generate_one_target_outcome_feature(outcome, raw_locations, edit_nucleotide, editing_window)
        outcomes.append(outcome)
        two_key_dic(outcome_features, seq, outcome, outcome_feature)
    return outcomes

def obtain_outcomes(seqs, editor, editing_window, target_nucleotide, edit_nucleotide):
    all_outcomes, outcome_features = {}, {}
    for seq_id, seq in seqs.items():
        outcomes = generate_one_target_outcomes(seq, outcome_features, target_nucleotide, edit_nucleotide, editing_window)
        all_outcomes[seq] = outcomes
    return all_outcomes, outcome_features

def merge_result(seqs, all_outcomes, outcome_features):
    df = []
    for seq_id, seq in seqs.items():
        for outcome in outcome_features[seq]:
            df.append([seq_id, seq, outcome, outcome_features[seq][outcome]])
    df = pd.DataFrame(df, columns = ['seq_id', 'target', 'outcome', 'outcome_feature'])
    return df

def get_nt_features(targets, seq_length):
    NTS = ['Aa','TtUu','Gg','Cc']
    x = np.zeros((len(targets), seq_length, len(NTS)), dtype=np.int32)
    for i, seq in enumerate(targets):
        for p, s in enumerate(seq):
            if s in NTS[0]: x[i][p][0] = 1
            elif s in NTS[1]: x[i][p][1] = 1
            elif s in NTS[2]: x[i][p][2] = 1
            elif s in NTS[3]: x[i][p][3] = 1
            else: raise Exception('Nucleotide not in {A,a,T,t,U,u,G,g,C,c,}; not accepted: %s'%s)
    return x

def get_outcome_features(df):
    outcome_features = np.array([[int(i) for i in list(feature)] for feature in df['outcome_feature'] ])
    return outcome_features

def get_binding_energies(df, filename):
    df_energy = pd.read_csv(filename, sep = '\t', index_col = 'guideSeq')
    x = np.array([[df_energy.loc[seq[4:27]]['CRISPRoff_score']] for seq in df['target']])
    return x

def get_cas9(df, filename):
    df_energy = pd.read_csv(filename, index_col = '30mer')
    x = np.array([[df_energy.loc[seq]['CRISPRon']] for seq in df['target']])
    return x

def get_weights(df, weight):
    x = np.array([weight for seq in df['target']])
    return x

def create_feature_X(df, off_filename, on_filename, seq_length, weight):
    x_seq = get_nt_features(df['target'], seq_length)
    x_outcome = get_outcome_features(df)
    x_energy = get_binding_energies(df, off_filename)
    x_cas9 = get_cas9(df, on_filename)
    x_weight = get_weights(df, weight)
    X = [x_seq, x_outcome, x_energy, x_cas9, x_weight]
    return X

def predict_X(path_to_model, editor, X, df):
    test_preds = []
    DFS = glob.glob(f'{path_to_model}/{editor}/model*')
    if len(DFS) < 1:
        raise Exception(f'{path_to_model}/{editor}/model* contains no models')
    for path in glob.glob(f'{path_to_model}/{editor}/model*'):
        model = keras.models.load_model(path)
        test_pred = model.predict(X)
        test_preds.append(test_pred)
    df_pred = np.average(np.array(test_preds), axis = 0)
    df_pred = pd.DataFrame(df_pred, columns = ['pred_eff', 'pred_freq'])
    df = pd.concat([df, df_pred], axis = 1)
    return df

def normalize_dataset(df):
    df['pred_freq'] = df['pred_freq'].apply(lambda x: x if x > 0 else 0)
    predicted_editingeffs = df.groupby(by=['seq_id']).mean(numeric_only = True)['pred_eff'].to_dict()
    df['pred_eff'] = df['seq_id'].apply(lambda x: predicted_editingeffs[x])
    predicted_sum_outcomefreqs = df.groupby(by=['seq_id']).sum(numeric_only = True)['pred_freq'].to_dict()
    df['pred_freq'] = df.apply(lambda row:
        row['pred_freq']/predicted_sum_outcomefreqs[row['seq_id']]*predicted_editingeffs[row['seq_id']], axis = 1)
    return df

def save_prediction_result(df, editor, path_to_out):
    output = f'{path_to_out}/crispron{editor}_prediction.tsv'
    df.drop(columns = ['outcome_feature'], inplace = True)
    df = normalize_dataset(df)
    df.sort_values(by=['seq_id', 'pred_freq'], inplace = True, ascending = False)
    #numeric_columns = df.select_dtypes(include='number').columns
    #df[numeric_columns] = df[numeric_columns].round(4)
    df.to_csv(output, sep = '\t', header = True, index = False)

def main():
    editor, fasta_filename, off_filename, on_filename, path_to_model, path_to_out = sys.argv[1:7]
    seq_length, editing_window, target_nucleotide, edit_nucleotide, weight = get_setting(editor)
    w = running_log(path_to_out)
    seqs = obtain_seqs(editor, fasta_filename, seq_length, editing_window, w)
    all_outcomes, outcome_features = obtain_outcomes(seqs, editor, editing_window, target_nucleotide, edit_nucleotide)
    df = merge_result(seqs, all_outcomes, outcome_features)
    X = create_feature_X(df, off_filename, on_filename, seq_length, weight)
    df = predict_X(path_to_model, editor, X, df)
    save_prediction_result(df, editor, path_to_out)
    w.close()

if __name__ == '__main__':
    main()
