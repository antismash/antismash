#!/usr/bin/env python3
# License: GNU Affero General Public License v3 or later
# A copy of GNU AGPL v3 should have been included in this software package in LICENSE.txt.

import importlib.util
import os

if __name__ == "__main__":
    path = os.path.join(__file__.split(os.path.join("antismash", "modules"))[0],
                        os.path.join("antismash", "common", "external", "rodeo_svm", "svm_classify.py"))

    spec = importlib.util.spec_from_file_location("svm_classify", path)
    svm_classify = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(svm_classify)

    data_dir = os.path.dirname(__file__)
    prefix = os.path.join(data_dir, "lanthipeptide")
    training_set = os.path.join(data_dir, "training_set.csv")

    svm_classify.save_classifier(training_set, prefix, kernel='rbf', C=9.77e6,
                                 gamma=1.78e-9)

    assert os.path.exists(prefix + ".scaler.pkl")
    assert os.path.exists(prefix + ".classifier.pkl")
    print("scaler and classifier regenerated")
    print("  located in:", data_dir)
