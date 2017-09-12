
'''
====CAML v0.1====
cschwal's assisted machine learning suite for RODEO

'''
print(__doc__)

import svm_sample
import svm_optimize

invalidinp = True

def main():
    print('Indicate which stage of SVM classification to perform:\n 1...splitting data set into training and fitting\n 2...optimizing hyperparameters\n 0...exit\n')

    stageinput = int(input('indicate task: \n'))

    if (stageinput == 1) or (stageinput == 2) or (stageinput == 0):
        if stageinput == 1:
            print('selecting training set...\n')
            svm_sample.main()
        elif stageinput == 2:
            print('optimizing hyperparameters...\n')
            svm_optimize.main()
        elif stageinput == 0:
            exit()
    else:
        print('INPUT ERROR\n')

while invalidinp:
    main()
