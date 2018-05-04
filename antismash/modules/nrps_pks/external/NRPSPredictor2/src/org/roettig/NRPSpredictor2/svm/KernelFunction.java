package org.roettig.NRPSpredictor2.svm;

public interface KernelFunction
{
	double compute(FeatureVector v1, FeatureVector v2);
}
