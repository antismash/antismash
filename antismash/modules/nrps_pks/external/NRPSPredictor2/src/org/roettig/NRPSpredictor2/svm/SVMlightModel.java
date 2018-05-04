package org.roettig.NRPSpredictor2.svm;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.ArrayList;
import java.util.List;

public class SVMlightModel
{
	public SVMlightModel(InputStream in) throws FileNotFoundException, IOException
	{
		parseModel(in);
	}
	
	private static double readFrontDouble(String line)
	{
		StringBuffer buf = new StringBuffer("");
		char c;
		int i = 0;
		while((c=line.charAt(i))!=' ')
		{
			buf.append(c);
			i++;
		}
		return Double.parseDouble(buf.toString());
	}
	
	private static int readFrontInt(String line)
	{
		StringBuffer buf = new StringBuffer("");
		char c;
		int i = 0;
		while((c=line.charAt(i))!=' ')
		{
			buf.append(c);
			i++;
		}
		return Integer.parseInt(buf.toString());
	}
	
	private void parseModel(InputStream in) throws IOException
	{
		BufferedReader reader = new BufferedReader(new InputStreamReader(in));
		String line;
		
		// read header
		line = reader.readLine();
		// read kernel type
		line = reader.readLine();
		kerneltype = readFrontInt(line);
		// read kernel param 1
		line = reader.readLine();
		kernelparam1 = readFrontDouble(line);
		// read kernel param 2
		line = reader.readLine();
		kernelparam2 = readFrontDouble(line);
		
		// ignore these lines
		line = reader.readLine();
		line = reader.readLine();
		line = reader.readLine();

		// read dim
		line = reader.readLine();
		dim  = readFrontInt(line);
		
		// ignore these lines
		line = reader.readLine();
		line = reader.readLine();
		
		// read bias
		line = reader.readLine();
		bias = readFrontDouble(line);


		while((line=reader.readLine())!=null)
		{
			// split "fvec # comment" by #
			String[] toks = line.split("#");
			SupportVector sv = readSV(toks[0],dim);
			svs.add(sv);
		}
		
	}
	
	private static SupportVector readSV(String line, int dim)
	{
		SupportVector sv = new SupportVector(dim);
		
		String[] toks = line.split("\\s+");
		double yalpha = Double.parseDouble(toks[0]);
		
		sv.setYalpha(yalpha);
		
		for(int i=1;i<toks.length;i++)
		{
			Double[]  val  = new Double[1];
			Integer[] idx  = new Integer[1];
			readFeature(toks[i],idx,val);
			sv.setFeatureValue(idx[0]-1, val[0]);
		}
		
		return sv;
	}
	
	private static void readFeature(String s, Integer[] idx, Double[] val)
	{
		String[] toks = s.split(":");
		idx[0] = Integer.parseInt(toks[0]);
		val[0] = Double.parseDouble(toks[1]);
	}
	
	public double predict(FeatureVector x, KernelFunction k_fun)
	{
		double ret = 0.0;
		for(SupportVector sv: svs)
		{
			ret += sv.getYalpha() * k_fun.compute(sv, x);
		}
		return ret-bias;
	}

	public double predict(FeatureVector x)
	{
		KernelFunction k_fun = new LinearKernel();
		if(this.kerneltype==2)
		{
			k_fun = new RBFKernel(kernelparam2);
		}
		
		double ret = 0.0;
		for(SupportVector sv: svs)
		{
			ret += sv.getYalpha() * k_fun.compute(sv, x);
		}
		return ret-bias;
	}
	
	
	public static List<FeatureVector> readFromFile(String filename) throws IOException
	{
		List<FeatureVector> fvecs = new ArrayList<FeatureVector>();
		
		BufferedReader reader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(filename))));
		String line;
		
		while((line=reader.readLine())!=null)
		{
			String[] toks = line.split("#");
			SupportVector sv = readSV(toks[0],408);
			fvecs.add(sv);
		}
		
		return fvecs;
	}

	private int     dim;
	private int     kerneltype;
	private double  bias;
	private double  kernelparam1;
	private double  kernelparam2;
	
	private List<SupportVector> svs = new ArrayList<SupportVector>();
	
	
	public static void main(String[] args) throws FileNotFoundException, IOException
	{
		SVMlightModel mdl = new SVMlightModel(new FileInputStream(new File("/tmp/model")));
		/*
		FeatureVector fv = new FeatureVector(408);
		fv.setFeatureValue(0, 1.0);
		fv.setFeatureValue(2, 1.0);
		fv.setFeatureValue(54, 1.0);
		System.out.println("y="+mdl.predict(fv, new RBFKernel(0.0056731971)));
		
		fv = new FeatureVector(408);
		fv.setFeatureValue(0, .0);
		fv.setFeatureValue(2, 0.0);
		fv.setFeatureValue(32, 0.0);
		System.out.println("y="+mdl.predict(fv, new RBFKernel(0.0056731971)));
		*/
		List<FeatureVector> fvecs = SVMlightModel.readFromFile("/tmp/test");
		for(FeatureVector fvec: fvecs)
		{
			//System.out.println("y="+mdl.predict(fvec, new RBFKernel(0.0056731971)));
			System.out.println("y="+mdl.predict(fvec));
		}
	}

}
