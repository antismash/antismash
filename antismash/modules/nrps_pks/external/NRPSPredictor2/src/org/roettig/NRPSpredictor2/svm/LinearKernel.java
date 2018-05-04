package org.roettig.NRPSpredictor2.svm;

public class LinearKernel implements KernelFunction
{

	@Override
	public double compute(FeatureVector v1, FeatureVector v2)
	{
		double[] f1 = v1.getFeatures();
		double[] f2 = v2.getFeatures();
		double sum = 0;
		for(int i=0;i<f1.length;i++)
			sum += f1[i]*f2[i];
		return sum;
	}

}
