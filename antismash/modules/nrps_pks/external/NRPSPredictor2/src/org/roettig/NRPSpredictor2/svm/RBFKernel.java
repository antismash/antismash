package org.roettig.NRPSpredictor2.svm;

public class RBFKernel implements KernelFunction
{

	private double gamma;
	
	public RBFKernel(double gamma)
	{
		this.gamma = gamma;
	}
	
	@Override
	public double compute(FeatureVector v1, FeatureVector v2)
	{
		double d = FeatureVector.dist(v1, v2);
		return Math.exp( -gamma*d*d);
	}
	
	public double getGamma()
	{
		return gamma;
	}

	public void setGamma(double gamma)
	{
		this.gamma = gamma;
	}

	
}
