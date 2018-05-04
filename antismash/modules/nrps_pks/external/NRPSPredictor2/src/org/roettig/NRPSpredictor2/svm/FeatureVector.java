package org.roettig.NRPSpredictor2.svm;

public class FeatureVector
{
	private double[] feats;
	
	public FeatureVector(double[] feats)
	{
		this.feats = feats;
	}
	
	public FeatureVector(int n)
	{
		this.feats = new double[n];
	}
	
	public double[] getFeatures()
	{
		return feats;
	}
	
	public void setFeatureValue(int idx, double value)
	{
		feats[idx] = value;
	}
	
	public double getFeatureValue(int idx)
	{
		return feats[idx];
	}
	
	public final static double dist(FeatureVector v1, FeatureVector v2)
	{
		double[] f1 = v1.getFeatures();
		double[] f2 = v2.getFeatures();
		double sum = 0.0;
		for(int i=0;i<f1.length;i++)
		{
			double d = f1[i]-f2[i]; 
			sum += d*d;
		}	
		return Math.sqrt(sum);
	}
}
