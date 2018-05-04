package org.roettig.NRPSpredictor2;

import java.beans.XMLDecoder;
import java.beans.XMLEncoder;
import java.io.BufferedInputStream;
import java.io.BufferedOutputStream;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;
import java.util.Locale;

import org.roettig.utilities.StringHelper;

public class ADomain implements Serializable
{
	private static final long	serialVersionUID	= -8665277739063399904L;
	
	public String sid;
	public String sig8a;
	public String sigstach;
	public String spec;
	public int startPos;
	public int endPos;
	public double pfamscore;
	public boolean outlier;
	
	public static String NRPS2_THREE_CLUSTER  = "NRPS2_THREE_CLUSTER";
	public static String NRPS2_THREE_CLUSTER_FUNGAL  = "NRPS2_THREE_CLUSTER_FUNGAL";
	public static String NRPS2_LARGE_CLUSTER  = "NRPS2_LARGE_CLUSTER";
	public static String NRPS2_SMALL_CLUSTER  = "NRPS2_SMALL_CLUSTER";
	public static String NRPS2_SINGLE_CLUSTER = "NRPS2_SINGLE_CLUSTER";
	public static String NRPS2_STACH_NN       = "NRPS2_STACH_NN";
	public static String NRPS2_STACH_NN_FUNGAL= "NRPS2_STACH_NN_FUNGAL";
	public static String NRPS1_LARGE_CLUSTER  = "NRPS1_LARGE_CLUSTER";
	public static String NRPS1_SMALL_CLUSTER  = "NRPS1_SMALL_CLUSTER";
	
	private List<Detection> detections = new ArrayList<Detection>();
	
	public boolean hasDetections(String pattern)
	{
		for(Detection det: detections)
		{
			if(det.getType().contains(pattern))
				return true;
		}
		return false;
	}
	
	public Detection getBestDetection(String type)
	{
		double maxScore = Double.NEGATIVE_INFINITY;
		Detection best  = null;
		for(Detection det: detections)
		{
			if(det.getType().equals(type))
			{
				if(det.getScore()>maxScore)
				{
					maxScore = det.getScore();
					best     = det;
				}
			}
		}
		return best;
	}
	
	public List<Detection> getDetections()
	{
		return detections;
	}

	public void setDetections(List<Detection> detections)
	{
		this.detections = detections;
	}

	public void addDetection(String type, String label, double score)
	{
		detections.add( new Detection(type,label,score) );
	}
	
	public void addDetection(String type, String label, double score, double prec)
	{
		detections.add( new Detection(type,label,score, prec) );
	}

	public String getSid()
	{
		return sid;
	}

	public void setSid(String sid)
	{
		this.sid = sid;
	}

	public String getSig8a()
	{
		return sig8a;
	}

	public void setSig8a(String sig8a)
	{
		this.sig8a = sig8a;
		if(sigstach==null)
		{
			StringBuffer sb = new StringBuffer();
			sb.append(sig8a.charAt(5));
			sb.append(sig8a.charAt(6));
			sb.append(sig8a.charAt(9));
			sb.append(sig8a.charAt(12));
			sb.append(sig8a.charAt(14));
			sb.append(sig8a.charAt(16));
			sb.append(sig8a.charAt(21));
			sb.append(sig8a.charAt(29));
			sb.append(sig8a.charAt(30));
			sb.append('K');
			sigstach = sb.toString();
		}
	}

	public String getSigstach()
	{
		return sigstach;
	}

	public void setSigstach(String sigstach)
	{
		this.sigstach = sigstach;
	}

	public String toString()
	{
		StringBuffer sb = new StringBuffer();
		sb.append(sid+"|");
		sb.append(sigstach+"|");
		sb.append(sig8a+"|");
		return sb.toString();
	}
	
	public void write(String filename) throws Exception
	{
        XMLEncoder encoder =
           new XMLEncoder(
              new BufferedOutputStream(
                new FileOutputStream(filename)));
        encoder.writeObject(this);
        encoder.close();
    }
	
	public static ADomain read(String filename) throws Exception 
	{
        XMLDecoder decoder =
            new XMLDecoder(new BufferedInputStream(
                new FileInputStream(filename)));
        ADomain o = (ADomain)decoder.readObject();
        decoder.close();
        return o;
    }

	public boolean isOutlier()
	{
		return outlier;
	}

	public void setOutlier(boolean outlier)
	{
		this.outlier = outlier;
	}
	
	
	
	public int getStartPos()
	{
		return startPos;
	}

	public void setStartPos(int startPos)
	{
		this.startPos = startPos;
	}

	public int getEndPos()
	{
		return endPos;
	}

	public void setEndPos(int endPos)
	{
		this.endPos = endPos;
	}

	public double getPfamscore()
	{
		return pfamscore;
	}

	public void setPfamscore(double pfamscore)
	{
		this.pfamscore = pfamscore;
	}

	public static void main(String[] args)
	{
		ADomain ad = new ADomain();
		ad.addDetection(ADomain.NRPS2_SINGLE_CLUSTER, "asp", 1.2);
		ad.addDetection(ADomain.NRPS2_SINGLE_CLUSTER, "glu", 1.5);
		ad.addDetection(ADomain.NRPS2_SINGLE_CLUSTER, "ala", 0.5);
		System.out.println(ad.getBestDetection(ADomain.NRPS2_SINGLE_CLUSTER));
	}


}