package org.roettig.NRPSpredictor2;

import gnu.getopt.Getopt;

import java.beans.XMLDecoder;
import java.beans.XMLEncoder;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.roettig.NRPSpredictor2.svm.FeatureVector;
import org.roettig.NRPSpredictor2.svm.SVMlightModel;

import libsvm.svm;
import libsvm.svm_model;
import libsvm.svm_node;

public class NRPSpredictor2
{	
	
	public static void main(String[] argv) throws Exception
	{
		try
		{
			datadir = System.getProperty("datadir",datadir);

			if(argv.length==0)
			{
				System.out.println("Usage: NRPSPredictor2 -i <inputfile> -r <reportfile> -s [0|1 use signatures?]\n");
				System.exit(0);
			}

			
			initSigDB();

			parseCommandline(argv);

			banner();

			if(useNRPS1input)
			{
				System.out.println("## using NRPS1 input");
				parseNRPS1(inputfile);
			}
			else
			{
				System.out.println("## using signature input");
				parseSigs(inputfile);
			}
			

			System.out.println("## start predicting on "+adoms.size()+" signatures");

			// we now have a list of adomain objects
			if(bacterialMode)
			{
				System.out.println("## bacterial mode");
				bacterialPrediction();
			}
			else
			{
				System.out.println("## fungal mode");
				fungalPrediction();
			}

			// export annotated adomains

			if(outputfile!=null)
				store(outputfile,adoms);
			if(reportfile!=null)
				report(reportfile,adoms);
		}
		catch(Throwable e)
		{
			e.printStackTrace();
			System.exit(1);
		}
	}
	
	public static void banner()
	{
		System.out.println("\n");
		System.out.println(" ##        Welcome to NRPSpredictor2 by        ##");
		System.out.println(" ##    Marc Roettig, Marnix Medema, Kai Blin   ##");
		System.out.println(" ##     based on work by Christian Rausch      ##\n");
		System.out.println(" please cite: http://dx.doi.org/10.1093/nar/gki885");
		System.out.println(" please cite: http://dx.doi.org/10.1093/nar/gkr323\n\n");
	}
	
	private static boolean useNRPS1input  = false;
	private static boolean bacterialMode  = true;
	private static String inputfile;
	private static String outputfile;
	private static String reportfile;
	private static String datadir;
	
	private static List<ADomain> adoms = new ArrayList<ADomain>();
	
	public static void parseSigs(String filename) throws Exception
	{
		if(!checkSignatureFormat(new File(filename)))
			crash("invalid signature file supplied");
		
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = null;
		ADomain cur_adom = null;
		while ((line = br.readLine()) != null)   
		{
			line = line.trim();
			if(line.equals(""))
				continue;
			String toks[] = line.split("\\\t");
			cur_adom = new ADomain();
			cur_adom.setSig8a(toks[0]);
			if(toks.length>1)
				cur_adom.sid = toks[1];
			adoms.add(cur_adom);
		}
	}
	
	public static boolean checkSignatureFormat(File infile) throws IOException
	{
		BufferedReader br = new BufferedReader(new FileReader(infile));
		String line = "";
		while((line=br.readLine())!=null)
		{
			line = line.trim();
			if(line.equals(""))
				continue;
			String toks[] = line.split("\\t");
			if(toks.length!=2)
			{
				return false;
			}
			String sig = toks[0].toUpperCase();
			if(sig.length()!=34)
			{
				return false;
			}
			if(!sig.matches("[A,C,D,E,F,G,H,I,K,L,M,N,P,Q,R,S,T,V,W,X,Y,-]+"))
				return false;
		}
		br.close();
		return true;
	}
	
	private static void crash(String msg)
	{
		System.err.println(msg);
		System.exit(1);
	}

	public static String fixLabel(String lab)
	{
		
		// large cluster predictions
		//String large_cluster[] = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"}; 
	
		// small cluster predictions
		//String small_cluster[] = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"};	 

		
		String ret = "N/A";
		if(lab.equals("phe=trp=phg=tyr=bht"))
			ret = "phe,trp,phg,tyr,bht";
		if(lab.equals("ser=thr=ser-thr=dht=dhpg=dpg=hpg"))
			ret = "ser,thr,dhpg,hpg";
		if(lab.equals("gly=ala=val=leu=ile=abu=iva"))
			ret = "gly,ala,val,leu,ile,abu,iva";
		if(lab.equals("asp=asn=glu=gln=aad"))
			ret = "asp,asn,glu,gln,aad";
		if(lab.equals("cys"))
			ret = "cys";
		if(lab.equals("orn=lys=arg"))
			ret = "orn,lys,arg";
		if(lab.equals("pro=pip"))
			ret = "pro,pip";
		if(lab.equals("dhb=sal"))
			ret = "dhb,sal";
		
		if(lab.equals("aaf"))
			ret = "aad";
		if(lab.equals("val=leu=ile=abu=iva"))
			ret = "val,leu,ile,abu,iva";
		if(lab.equals("arg"))
			ret = "arg";		
		if(lab.equals("asp=asn"))
			ret = "asp,asn";
		if(lab.equals("cys"))
			ret = "cys";
		if(lab.equals("dhb=sal"))
			ret = "dhb,sal";
		if(lab.equals("glu=gln"))
			ret = "glu,gln";
		if(lab.equals("orn"))
			ret = "orn,horn";		
		if(lab.equals("tyr=bht"))
			ret = "tyr,bht";		
		if(lab.equals("pro"))
			ret = "pro";		
		if(lab.equals("ser"))
			ret = "ser";		
		if(lab.equals("dhpg=dpg=hpg"))
			ret = "dhpg,hpg";		
		if(lab.equals("phe=trp"))
			ret = "phe,trp";
		if(lab.equals("gly=ala"))
			ret = "gly,ala";		
		if(lab.equals("thr=dht"))
			ret = "thr,dht";		
		return ret;
	}
	
	public static void parseNRPS1(String filename) throws Exception
	{
		BufferedReader br = new BufferedReader(new FileReader(filename));
		String line = null;
		ADomain cur_adom = null;
		while ((line = br.readLine()) != null)   
		{
			if(line.startsWith("[no hits above thresholds]"))
			{
				cur_adom = null;
			}
			if(line.startsWith("Score: "))
			{
				String scoreline = line.substring(6).trim();
				double score = 0;
				try
				{
					score = Double.parseDouble(scoreline);
				}
				catch(Exception e)
				{
					
				}
				cur_adom.pfamscore = score;
			}
			if(line.startsWith("Loc: "))
			{
				String locsline = line.substring(4).trim();
				String[] locs = locsline.split("_");
				
				int start = -1;
				int end   = -1;
				try
				{
					start = Integer.parseInt(locs[0]);
					end   = Integer.parseInt(locs[1]);
				}
				catch(Exception e)
				{
					
				}
				cur_adom.startPos = start;
				cur_adom.endPos   = end;
			}
			if(line.startsWith("Modul: "))
			{
				if(cur_adom!=null)
					adoms.add(cur_adom);
				cur_adom = new ADomain();

				String sid = line.substring(7).trim();
				cur_adom.sid = sid;
				continue;
			}
			if(line.contains("10 amino acid code"))
			{
				String toks[] = line.split(":");
				String sig = toks[1].trim();
				cur_adom.sigstach = sig;
				continue;
			}
			if(line.contains("8A-Code"))
			{
				String toks[] = line.split(":");
				String sig = toks[1].trim();
				cur_adom.sig8a = sig;
				continue;
			}
			if(line.contains("large clusters"))
			{
				while((line = br.readLine())!=null)
				{
					if(line.contains("small clusters"))
						break;
					if(line.contains("[no predictions"))
					{
						continue;
					}
					String toks[] = line.split(":");
					String label = fixLabel(toks[0]);
					if(precsNRPS1.containsKey(label))
						cur_adom.addDetection(ADomain.NRPS1_LARGE_CLUSTER, toks[0], Double.parseDouble(toks[1]),precsNRPS1.get(label));
					else
						cur_adom.addDetection(ADomain.NRPS1_LARGE_CLUSTER, toks[0], Double.parseDouble(toks[1]),0.0);
				}
				while((line = br.readLine())!=null)
				{
					if(line.contains("Alis (>="))
						break;
					if(line.contains("[no predictions"))
					{
						continue;
					}
					String toks[] = line.split(":");
					String label = fixLabel(toks[0]);
					if(precsNRPS1.containsKey(label))
						cur_adom.addDetection(ADomain.NRPS1_SMALL_CLUSTER, toks[0], Double.parseDouble(toks[1]),precsNRPS1.get(label));
					else
						cur_adom.addDetection(ADomain.NRPS1_SMALL_CLUSTER, toks[0], Double.parseDouble(toks[1]),0.0);
					//cur_adom.addDetection(ADomain.NRPS1_SMALL_CLUSTER, toks[0], Double.parseDouble(toks[1]));
				}
				continue;
			}
			
		}
		
		if(cur_adom!=null)
			adoms.add(cur_adom);
		
		br.close();
		
		//store(outputdir+"/preds",adoms);
		
	}
	
	public static void store(String filename, List<ADomain> data) throws IOException
	{
		XMLEncoder encoder =
	           new XMLEncoder(
	              new BufferedOutputStream(
	                new FileOutputStream(outputfile)));
	    encoder.writeObject(data);
	    encoder.close();
	}
	
	public static List<ADomain> load(String filename) throws Exception
	{
		XMLDecoder decoder =
            new XMLDecoder(new BufferedInputStream(
                new FileInputStream(filename)));
        List<ADomain> o = (List<ADomain>) decoder.readObject();
        decoder.close();
        return o;
	}
	
	
	
	private static List<String> specs     = new ArrayList<String>();
	private static List<String> specs_fungal = new ArrayList<String>();
	private static List<String> sig8a     = new ArrayList<String>();
	private static List<String> sigstach = new ArrayList<String>();
	private static List<String> sigstach_fungal = new ArrayList<String>();
	
	private static void initSigDB() throws IOException
	{
		BufferedReader in = new BufferedReader(new FileReader(datadir+"/labeled_sigs"));
		String line = "";
		while((line=in.readLine())!=null)
		{
			String toks[] = line.split("\\\t");
			specs.add(toks[0]);
			sig8a.add(toks[1]);
			sigstach.add(toks[2]);
		}
		in.close();
		
		in = new BufferedReader(new FileReader(datadir+"/fungal_labeled_sigs"));
		line = "";
		while((line=in.readLine())!=null)
		{
			String toks[] = line.split("\\\t");
			specs_fungal.add(toks[0]);
			sigstach_fungal.add(toks[2]);
		}
		in.close();
	}
	
	

	private static void report(String outputfile, List<ADomain> adoms) throws FileNotFoundException
	{
		PrintWriter out = new PrintWriter( new File(reportfile) );
		out.println(String.format("#%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s<tab>%s","sequence-id","8A-signature","stachelhaus-code","3class-pred","large-class-pred","small-class-pred","single-class-pred","nearest stachelhaus code","NRPS1pred-large-class-pred","NRPS2pred-large-class-pred","outside applicability domain?","coords","pfam-score"));
		for(ADomain adom: adoms)
		{
			// id <tab> 8a_sig <tab> stachelhaus_sig <tab> three_class_prediction_nrps2 <tab> large_cluster_prediction_nrps2 <tab> small_cluster_prediction_nrps2 <tab> single_aa_prediction_nrps2 <tab> nearest_stachelhaus_specificity <tab> large_cluster_prediction_nrps1 <tab> small_cluster_prediction_nrps1 <tab> isAdomainUnusual? <tab> coords <pfam-score>
			
			Detection best_THREE_CLASS = adom.getBestDetection(ADomain.NRPS2_THREE_CLUSTER);
			String pred_THREE_CLASS = (best_THREE_CLASS!=null?best_THREE_CLASS.getLabel():"N/A");
			
			Detection best_LARGE_CLASS = adom.getBestDetection(ADomain.NRPS2_LARGE_CLUSTER);
			String pred_LARGE_CLASS = (best_LARGE_CLASS!=null?best_LARGE_CLASS.getLabel():"N/A");
			
			Detection best_SMALL_CLASS = adom.getBestDetection(ADomain.NRPS2_SMALL_CLUSTER);
			String pred_SMALL_CLASS = (best_SMALL_CLASS!=null?best_SMALL_CLASS.getLabel():"N/A");
			
			Detection best_SINGLE_CLASS = adom.getBestDetection(ADomain.NRPS2_SINGLE_CLUSTER);
			String pred_SINGLE_CLASS = (best_SINGLE_CLASS!=null?best_SINGLE_CLASS.getLabel():"N/A");
			
			Detection best_STACH_NN = adom.getBestDetection(ADomain.NRPS2_STACH_NN);
			String pred_STACH_NN = (best_STACH_NN!=null?best_STACH_NN.getLabel():"N/A");
			
			Detection best_LARGE_CLASS1 = adom.getBestDetection(ADomain.NRPS1_LARGE_CLUSTER);
			String pred_LARGE_CLASS1 = (best_LARGE_CLASS1!=null?best_LARGE_CLASS1.getLabel():"N/A");
			
			Detection best_SMALL_CLASS1 = adom.getBestDetection(ADomain.NRPS1_SMALL_CLUSTER);
			String pred_SMALL_CLASS1 = (best_SMALL_CLASS1!=null?best_SMALL_CLASS1.getLabel():"N/A");
			
			out.println(String.format("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%e",adom.sid,adom.sig8a,adom.sigstach, pred_THREE_CLASS, pred_LARGE_CLASS, pred_SMALL_CLASS, pred_SINGLE_CLASS, pred_STACH_NN, pred_LARGE_CLASS1, pred_SMALL_CLASS1 ,(adom.isOutlier()?1:0),String.format("%d:%d",adom.startPos,adom.endPos),adom.pfamscore));
		}
		out.close();
	}
	
	public static void fungalPrediction() throws ParseException, IOException
	{	
		PrimalEncoder encR = new PrimalRauschEncoder();
		
		// check applicability domain
		//checkAD("bacterial");
		
		//////////
		// NRPS1
		/////////
		
		// large cluster predictions
		String large_cluster[] = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"}; 
	
		// small cluster predictions
		String small_cluster[] = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"};	 
	
		
		if(!useNRPS1input)
		{
			for(String sc: large_cluster)
			{
				detect(ADomain.NRPS1_LARGE_CLUSTER, sc, encR);
			}

			for(String sc: small_cluster)
			{			
				detect(ADomain.NRPS1_SMALL_CLUSTER, sc, encR);
			}
		}
		
		//////////
		// NRPS2 
		/////////
	
		
		PrimalEncoder encW = new PrimalWoldEncoder();
		
		// 3 class predictions
		String three_class[] = {"hydrophilic","hydrophobic-aliphatic","hydrophobic-aromatic"};
		for(String sc: three_class)
		{
			detect(ADomain.NRPS2_THREE_CLUSTER, sc, encW);
		}
		
		// large cluster predictions
		//String large_cluster[] = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"}; 
		
		for(String sc: large_cluster)
		{
			detect(ADomain.NRPS2_LARGE_CLUSTER, sc, encW);
		}
		
		// small cluster predictions
		//String small_cluster[] = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"}; 
		for(String sc: small_cluster)
		{			
			detect(ADomain.NRPS2_SMALL_CLUSTER, sc, encW);
		}
		
		// single aa predictions
		String single_cluster[] = {"aad","ala","arg","asn","asp","bht","cys","dhb","dhpg","gln","glu","gly","hpg","ile","iva","leu","lys","orn","phe","pip","pro","ser","thr","trp","tyr","val"}; 
		for(String sc: single_cluster)
		{
			detect(ADomain.NRPS2_SINGLE_CLUSTER, sc, encW);
		}
		
		detectStachNN(false);
	}
	
	public static void bacterialPrediction() throws ParseException, IOException
	{	
		PrimalEncoder encR = new PrimalRauschEncoder();
		
		// check applicability domain
		checkAD("bacterial");
		
		//////////
		// NRPS1
		/////////
		
		// large cluster predictions
		String large_cluster[] = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"}; 
	
		// small cluster predictions
		String small_cluster[] = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"};	 
	
		
		if(!useNRPS1input)
		{
			for(String sc: large_cluster)
			{
				detect(ADomain.NRPS1_LARGE_CLUSTER, sc, encR);
			}

			for(String sc: small_cluster)
			{			
				detect(ADomain.NRPS1_SMALL_CLUSTER, sc, encR);
			}
		}
		
		//////////
		// NRPS2 
		/////////
	
		
		PrimalEncoder encW = new PrimalWoldEncoder();
		
		// 3 class predictions
		String three_class[] = {"hydrophilic","hydrophobic-aliphatic","hydrophobic-aromatic"};
		for(String sc: three_class)
		{
			detect(ADomain.NRPS2_THREE_CLUSTER, sc, encW);
		}
		
		// large cluster predictions
		//String large_cluster[] = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"}; 
		
		for(String sc: large_cluster)
		{
			detect(ADomain.NRPS2_LARGE_CLUSTER, sc, encW);
		}
		
		// small cluster predictions
		//String small_cluster[] = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"}; 
		for(String sc: small_cluster)
		{			
			detect(ADomain.NRPS2_SMALL_CLUSTER, sc, encW);
		}
		
		// single aa predictions
		String single_cluster[] = {"aad","ala","arg","asn","asp","bht","cys","dhb","dhpg","gln","glu","gly","hpg","ile","iva","leu","lys","orn","phe","pip","pro","ser","thr","trp","tyr","val"}; 
		for(String sc: single_cluster)
		{
			detect(ADomain.NRPS2_SINGLE_CLUSTER, sc, encW);
		}
		
		detectStachNN(true);
	}
	
	public static void checkAD(String model) throws IOException
	{
		svm_model m = svm.svm_load_model(datadir+String.format("/models/%s_1class.mdl",model));
		
		PrimalWoldEncoder enc = new PrimalWoldEncoder();
		
		for(ADomain ad: adoms)
		{
			double fv[]  = enc.encode(ad.sig8a);
			svm_node x[] = makeSVMnode(fv,1);
			double yp = svm.svm_predict(m,x);
			if(yp>=0.0)
				ad.setOutlier(false);
			else
				ad.setOutlier(true);
		}
		
	}
	
	public static Map<String,Double> precs      = new HashMap<String,Double>();
	public static Map<String,Double> precsNRPS1 = new HashMap<String,Double>();
	
	public static void fillPrecs()
	{
		String three_class[]    = {"hydrophilic","hydrophobic-aliphatic","hydrophobic-aromatic"};
		String large_cluster[]  = {"phe,trp,phg,tyr,bht","ser,thr,dhpg,hpg","gly,ala,val,leu,ile,abu,iva","asp,asn,glu,gln,aad","cys","orn,lys,arg","pro,pip","dhb,sal"};
		String small_cluster[]  = {"aad","val,leu,ile,abu,iva","arg","asp,asn","cys","dhb,sal","glu,gln","orn,horn","tyr,bht","pro","ser","dhpg,hpg","phe,trp","gly,ala","thr"};
		String single_cluster[] = {"aad","ala","arg","asn","asp","bht","cys","dhb","dhpg","gln","glu","gly","hpg","ile","iva","leu","lys","orn","phe","pip","pro","ser","thr","trp","tyr","val"};
		
		precs.put("hydrophilic",0.940);
		precs.put("hydrophobic-aliphatic",0.974);
		precs.put("hydrophobic-aromatic",0.890);
		
		precs.put("phe,trp,phg,tyr,bht", 0.881);
		precs.put("ser,thr,dhpg,hpg", 0.967 );
		precs.put("gly,ala,val,leu,ile,abu,iva",0.947 );
		precs.put("asp,asn,glu,gln,aad", 0.969);
		precs.put("cys", 0.975);
		precs.put("orn,lys,arg", 0.898);
		precs.put("pro,pip", 0.867);
		precs.put("dhb,sal", 1.00);
		
		precs.put("aad",1.00 );
		precs.put("val,leu,ile,abu,iva", 0.892);
		precs.put("arg",1.00 );
		precs.put("asp,asn", 0.969);
		precs.put("cys", 0.983);
		precs.put("dhb,sal", 1.000);
		precs.put("glu,gln", 0.85);
		precs.put("orn,horn",0.900 );
		precs.put("tyr,bht", 0.892);
		precs.put("pro", 0.938);
		precs.put("ser", 1.0);
		precs.put("dhpg,hpg", 1.000);
		precs.put("phe,trp", 0.608);
		precs.put("gly,ala", 0.938);
		precs.put("thr", 0.978);
		
		precs.put("ala",0.901);
		precs.put("asn",0.934);
		precs.put("asp",0.700);
		precs.put("bht",0.782);
		precs.put("dhb",1.00);
		precs.put("dhpg",0.967);
		precs.put("gln",0.775);
		precs.put("glu",0.760);
		precs.put("gly",0.902);
		precs.put("hpg",1.0);
		precs.put("ile",1.0);
		precs.put("iva",0.933);
		precs.put("leu",0.957);
		precs.put("lys",0.500);
		precs.put("phe",0.740);
		precs.put("trp",0.400);
		precs.put("tyr",0.671);
		precs.put("val",0.801);
		
		precsNRPS1.put("phe,trp,phg,tyr,bht", 0.881);
		precsNRPS1.put("ser,thr,dhpg,hpg", 0.963 );
		precsNRPS1.put("gly,ala,val,leu,ile,abu,iva",0.940);
		precsNRPS1.put("asp,asn,glu,gln,aad", 0.969);
		precsNRPS1.put("cys", 0.958);
		precsNRPS1.put("orn,lys,arg", 0.898);
		precsNRPS1.put("pro,pip", 0.811);
		precsNRPS1.put("dhb,sal", 1.00);
		
		precsNRPS1.put("aad",1.00 );
		precsNRPS1.put("val,leu,ile,abu,iva", 0.90);
		precsNRPS1.put("arg",1.00 );
		precsNRPS1.put("asp,asn", 0.969);
		precsNRPS1.put("cys", 1.00);
		precsNRPS1.put("dhb,sal", 1.000);
		precsNRPS1.put("glu,gln", 0.86);
		precsNRPS1.put("orn,horn",0.80 );
		precsNRPS1.put("tyr,bht", 0.825);
		precsNRPS1.put("pro", 0.900);
		precsNRPS1.put("ser", 0.936);
		precsNRPS1.put("thr", 0.942);
		precsNRPS1.put("dhpg,hpg", 0.985);
		precsNRPS1.put("phe,trp", 0.671);
		precsNRPS1.put("gly,ala", 0.859);
		precsNRPS1.put("thr,dht", 0.942);
	}
	
	public static void detect(String type, String label, PrimalEncoder enc) throws ParseException, FileNotFoundException, IOException
	{
		SVMlightModel m = new SVMlightModel(new FileInputStream(String.format("%s/models/%s/[%s].mdl",datadir,type,label)));
		
		for(ADomain ad: adoms)
		{
			double fv[] = enc.encode(ad.sig8a);
			double yp   = m.predict(makeFVec(fv));
			if(yp>0.0)
			{
				if(type.contains("NRPS2"))
				{
					if(precs.containsKey(label))
						ad.addDetection(type, label, yp, precs.get(label));
					else
						ad.addDetection(type, label, yp, 0.0);
				}
				else
				{
					if(precsNRPS1.containsKey(label))
						ad.addDetection(type, label, yp, precsNRPS1.get(label));
					else
						ad.addDetection(type, label, yp, 0.0);
				}
			}
		}
	}
	
	public static void detectStachNN(boolean bacterial)
	{
		for(ADomain ad: adoms)
		{
			String[] hitspec = new String[1];
			double sim = nearestSignature(ad, hitspec, bacterial);
			if(bacterial)
				ad.addDetection(ADomain.NRPS2_STACH_NN, hitspec[0], sim);
			else
				ad.addDetection(ADomain.NRPS2_STACH_NN_FUNGAL, hitspec[0], sim);
		}
	}
	
	public static double nearestSignature(ADomain ad, String[] hit, boolean bacterial)
	{
		String a_sig   = ad.sigstach;
		int maxMatches = 0;
		int idx = 0;
		
		List<String> sigs;
		if(bacterial)
			sigs = sigstach;
		else
			sigs = sigstach_fungal;
		
		for(String sig: sigs)
		{
			int matches = 0;
			for(int i=0;i<sig.length();i++)
			{
				if(sig.charAt(i)==a_sig.charAt(i))
					matches++;
			}
			if(matches>maxMatches)
			{
				maxMatches = matches;
				hit[0] = specs.get(idx); 
			}
			idx++;
		}
		return (1.0*maxMatches)/(a_sig.length());
	}
	

	public static svm_node[] makeSVMnode(double[] row, int idx)
	{
		int m        = row.length;
		svm_node[] x = new svm_node[m+1];

		x[0]         = new svm_node();
		x[0].index   = 0; 
		x[0].value   = idx;

		for(int j=1;j<=m;j++)
		{
			x[j]       = new svm_node();
			x[j].index = j; 
			x[j].value = row[j-1];
		}
		return x;
	}
	
	public static FeatureVector makeFVec(double[] fts)
	{
		return new FeatureVector(fts);
	}
	
	public static void parseCommandline(String[] argv)
	{
		fillPrecs();
		
		Getopt g = new Getopt("NRPSpredictor2", argv, "i:o:b:s:r:");
		//
		int c;
		String arg;
		while ((c = g.getopt()) != -1)
		{
			switch(c)
			{
			case 's':
				arg = g.getOptarg();
				if(arg.equals("1"))
					useNRPS1input = false;
				else
					useNRPS1input = true;
				break;
				//
			case 'b':
				arg = g.getOptarg();
				if(arg.equals("1"))
					bacterialMode = true;
				else
					bacterialMode = false;
				break;
				//				
			case 'i':
				arg = g.getOptarg();
				inputfile = arg;
				break;
				//
			case 'r':
				arg = g.getOptarg();
				reportfile = arg;
				break;
				//				
			case 'o':
				arg = g.getOptarg();
				outputfile = arg;
				break;
				//
			default:
				break;				
			}
		}		
	}
	

}
