# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

""" Functions for rebuilding a SVM classifier """

import logging
import os

from .svm_classify import save_classifier


def pickle_classifier(training_set: str, prefix: str, overwrite: bool = False,
                      kernel: str = "rbf", C: float = 9.77e6, gamma: float = 1.78e-9) -> bool:
    """ Builds and saves a classifier in a pickled format, in the same directory as the training set.
        Files that are created:
            PREFIX.classifier.pkl
            PREFIX.scaler.pkl

        Args:
            training_set: the CSV file containing the training set
            prefix: a string to use as a prefix for the pickled files
            overwrite: if False and existing pickled classifiers exist, new ones will not be generated
            kernel: the kernel to use for scikits classifier
            C: the penalty parameter to use for the classifier
            gamma: the kernel coefficient for the classifier

        Returns:
            True if classifiers were written/updated, otherwise False
    """
    prefix = os.path.join(os.path.dirname(training_set), prefix)

    if not overwrite and os.path.exists(prefix + ".scaler.pkl") and os.path.exists(prefix + ".classifier.pkl"):
        return False

    logging.debug("rebuilding RODEO classifiers")
    save_classifier(training_set, prefix, kernel=kernel, C=C, gamma=gamma)

    if not os.path.exists(prefix + ".scaler.pkl") or not os.path.exists(prefix + ".classifier.pkl"):
        raise ValueError("rebuilding RODEO classifiers failed")
    return True
