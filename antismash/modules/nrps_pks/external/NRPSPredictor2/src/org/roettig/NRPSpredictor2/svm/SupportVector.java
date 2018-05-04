package org.roettig.NRPSpredictor2.svm;


public class SupportVector extends FeatureVector
{

	private double yalpha;
	
	public SupportVector(double[] feats, double yalpha)
	{
		super(feats);
		this.yalpha = yalpha;
	}
	
	public SupportVector(int n)
	{
		super(n);
	}
	
	public double getYalpha()
	{
		return yalpha;
	}
	
	public void setYalpha(double yalpha)
	{
		this.yalpha = yalpha;
	}
}
