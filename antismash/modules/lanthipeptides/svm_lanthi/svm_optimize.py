#! /usr/bin/env python
#
# Copyright (C) 2015 Jonathan I. Tietz
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Christopher J. Schwalen
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Parth S. Patel
# University of Illinois
# Department of Chemistry
#
# Copyright (C) 2015 Douglas A. Mitchell
# University of Illinois
# Department of Chemistry
#
# License: GNU Affero General Public License v3 or later
# Complete license availabel in the accompanying LICENSE.txt.
# or <http://www.gnu.org/licenses/>.
#
# This file is part of RODEO.
#
# RODEO is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# RODEO is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Affero General Public License for more details.
#
'''
   SVM optimization script

   Required input:
     a training file CSV

   Output:
     a CSV list of identifiers and classifications

   Note that all options must be hard-coded (the script does not take command-line arguments).

   RECOMMENDATION:
     Input CSV should ideally have its primary key as column 0, its classification as column 1, and feature as columns [2,...,end]
     For the fitting CSV, there will be no classification; leave it blank if you want (it'll be ignored upon import)
'''
import numpy as np
import csv

from sklearn import svm
from sklearn import preprocessing
from sklearn import cross_validation

# CONFIGURATION OPTIONS
# change these as desired

input_training_file = 'training_set.csv'         # the CSV containing the training set
output_filename = 'optimization_results.csv'     # output filename; this will be a CSV with the parameters and accuracy

primary_key_column = 0            # the column of the CSV that contains the primary key (identifier) for each record
classification_column = 1         # the column of the CSV that contains the classification for each record
csv_has_header = True             # set to true if CSVs have header rows to ignore upon import

# values for optimization
kernel_option = 'rbf'       # the options below are for rbf; feel free to modify for poly, linear, etc
C_max = 10
C_min = -1
C_base = 5
C_steps = 11
C_options = np.logspace(C_min,C_max,num=C_steps,base=C_base)
gamma_max = -2
gamma_min = -10
gamma_base = 10
gamma_steps = 9
gamma_options = np.logspace(gamma_min,gamma_max,num=gamma_steps,base=gamma_base)
class_weight_option = 'balanced'
folds_validation = [5]

def parse_CSV_to_dataset(csv_filename, dataset_type):
    '''Parse an input CSV into a data set

       Inputs:
            csv_filename            name of CSV file to be parsed
            dataset_type            either 'training' or 'fitting'
    '''
    dataset = []
    with open(csv_filename, 'rb') as csvfile:
      csv_read = csv.reader(csvfile, delimiter=',', quotechar='"')
      if csv_has_header:
        next(csv_read,None)

      if classification_column < primary_key_column:
        pk_first = 0
      else:
        pk_first = 1

      for row in csv_read:
        temp_entry = []
        temp_entry.append(row.pop(primary_key_column))
        temp_entry.append(row.pop(classification_column - pk_first))
        for c in row:
          temp_entry.append(float(c))
        # remove all unclassified features if training
        if dataset_type == 'training':
          if int(temp_entry[1]) == 1 or int(temp_entry[1]) == 0:
            temp_entry[1] = int(temp_entry[1])
            dataset.append(temp_entry)
        if dataset_type == 'fitting':
          dataset.append(temp_entry)
    return dataset

def write_to_csv(result_list, output_file):
    with open(output_filename, 'wb') as csvfile:
      csv_write = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
      csv_write.writerow(['kernel','fold cross-validation','C value','gamma value','class_weight','precision','recall','f1','score'])
      for row in result_list:
        csv_write.writerow(row)
    return

def main():

    # parse data

    print("Importing training and fitting data ...")
    training_data_unrefined = parse_CSV_to_dataset(input_training_file, 'training')

    training_data_just_features = []
    training_data_classifications = []

    for entry in training_data_unrefined:
        training_data_classifications.append(entry.pop(1))
        entry.pop(0)
        training_data_just_features.append(entry)
    print("  (done)")
    print("Initiating learning and fitting")

    # Scaling -- this ensures standardization of model and target data
    training_data_refined = preprocessing.scale(training_data_just_features)
    scaler = preprocessing.StandardScaler().fit(training_data_refined)
    training_data_refined = scaler.transform(training_data_just_features)

    test_results = []

    for fold in folds_validation:
        for C_option in C_options:
            for gamma_option in gamma_options:
                clf = svm.SVC(kernel=kernel_option,class_weight=class_weight_option,C=C_option,gamma=gamma_option)
                prec = cross_validation.cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='precision')
                recd = cross_validation.cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='recall')
                f1we = cross_validation.cross_val_score(clf, training_data_refined, training_data_classifications, cv=fold, scoring='f1')
                scor = prec.mean() * recd.mean()
                print('Using %d fold, %0.2E C, %0.2E gamma: %0.2f precision, %0.2f recall, %0.2f f1, %0.2f score' % (fold, C_option, gamma_option, prec.mean(), recd.mean(), f1we.mean(), scor))
                test_results.append([kernel_option,fold,C_option,gamma_option,class_weight_option,prec.mean(),recd.mean(),f1we.mean(),scor])
    maxscor = max(np.array(test_results)[:,-1])
    maxscor = float(maxscor)

    # Output results to file
    write_to_csv(test_results, output_filename)
    print(" ... Done")
    print('Maximum score of %0.5f achieved\n' % (maxscor))
    return

if __name__ == '__main__':
  main()
