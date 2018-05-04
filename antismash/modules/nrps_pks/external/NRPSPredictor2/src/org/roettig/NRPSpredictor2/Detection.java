package org.roettig.NRPSpredictor2;

import java.io.Serializable;
import java.util.Locale;

public class Detection  implements Serializable
{
	private static final long	serialVersionUID	= -5102112135793201058L;

	public Detection()
	{
		
	}
	
	public Detection(String _type, String _label, double _score)
	{
		type  = _type;
		label = _label;
		score = _score;
	}
	
	public Detection(String _type, String _label, double _score, double prec)
	{
		type  = _type;
		label = _label;
		score = _score;
		precision = prec;
	}
	
	private String type;
	private String label;
	private double precision;
	private double score;

	public String getType()
	{
		return type;
	}

	public void setType(String type)
	{
		this.type = type;
	}

	public String getLabel()
	{
		return label;
	}

	public void setLabel(String label)
	{
		this.label = label;
	}

	public double getScore()
	{
		return score;
	}

	public void setScore(double score)
	{
		this.score = score;
	}

	public String toString()
	{
		return String.format(Locale.ENGLISH,"[%s / %s] score: %.8f", type, label, score );
	}

	public double getPrecision()
	{
		return precision;
	}

	public void setPrecision(double precision)
	{
		this.precision = precision;
	}
	
	
}