
import numpy as np
import csv
# sklearn does things differently and pylint doesn't like that
# pylint: disable=no-name-in-module
from sklearn.model_selection import train_test_split
# pylint: enable=no-name-in-module

# CONFIGURATION OPTIONS
''' change these as desired '''

csv_filename = 'full_set.csv'         # the CSV containing the training set
dataset_type = 'full'
training_file = 'training_set.csv'
fitting_file = 'fitting_set.csv'

primary_key_column = 0;            # the column of the CSV that contains the primary key (identifier) for each record
classification_column = 1;         # the column of the CSV that contains the classification for each record
csv_has_header = True;             # set to true if CSVs have header rows to ignore upon import


def parse_CSV_to_dataset(csv_filename, dataset_type):
    dataset=[]
    with open(csv_filename, 'rb') as csvfile:
        csv_read = csv.reader(csvfile, delimiter=',', quotechar='"')
        for row in csv_read:
            dataset.append(row)
    if csv_has_header == True:
        header = dataset[0][:]
        dataset.pop(0)
    return dataset, header

def write_to_csv(result_list, output_file, header):
    with open(output_file, 'wb') as csvfile:
        csv_write = csv.writer(csvfile, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        if csv_has_header == True:
            result_list = np.vstack((header, result_list))
        for row in result_list:
            csv_write.writerow(row)
    return

def sampler(inpdata):
    datar = np.asarray(inpdata)
#    numfeats = datar.shape[1]-2
#    feats = datar[:,(numfeats-2):]
    classif = datar[:,classification_column]
    x_train, x_test, y_train, y_test = train_test_split(datar, classif, test_size=0.75, random_state=42)
    return x_train, x_test, y_train, y_test

def main():
    dataset, header=parse_CSV_to_dataset(csv_filename, dataset_type)
    x_train, x_test, y_train, y_test=sampler(dataset)
    write_to_csv(x_test, fitting_file, header)
    write_to_csv(x_train, training_file, header)

if __name__ == '__main__':
  main()
