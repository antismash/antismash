package org.roettig.NRPSpredictor2;


public class PrimalWoldEncoder implements PrimalEncoder
{

	private static double[] Z1 = new double[26];
	private static double[] Z2 = new double[26];
	private static double[] Z3 = new double[26];
	
	static
	{
		Z1[0]=0.07;Z1[17]=2.88;Z1[13]=3.22;Z1[3]=3.64;Z1[2]=0.71;Z1[16]=2.18;Z1[4]=3.08;Z1[6]=2.23;Z1[7]=2.41;Z1[8]=-4.44;Z1[11]=-4.19;Z1[10]=2.84;Z1[12]=-2.49;Z1[5]=-4.92;Z1[15]=-1.22;Z1[18]=1.96;Z1[19]=0.92;Z1[22]=-4.75;Z1[24]=-1.39;Z1[21]=-2.69;
		Z2[0]=-1.73;Z2[17]=2.52;Z2[13]=1.45;Z2[3]=1.13;Z2[2]=-0.97;Z2[16]=0.53;Z2[4]=0.39;Z2[6]=-5.36;Z2[7]=1.74;Z2[8]=-1.68;Z2[11]=-1.03;Z2[10]=1.41;Z2[12]=-0.27;Z2[5]=1.3;Z2[15]=0.88;Z2[18]=-1.63;Z2[19]=-2.09;Z2[22]=3.65;Z2[24]=2.32;Z2[21]=-2.53;
		Z3[0]=0.09;Z3[17]=-3.44;Z3[13]=0.84;Z3[3]=2.36;Z3[2]=4.13;Z3[16]=-1.14;Z3[4]=-0.07;Z3[6]=0.3;Z3[7]=1.11;Z3[8]=-1.03;Z3[11]=-0.98;Z3[10]=-3.14;Z3[12]=-0.41;Z3[5]=0.45;Z3[15]=2.23;Z3[18]=0.57;Z3[19]=-1.4;Z3[22]=0.85;Z3[24]=0.01;Z3[21]=-1.29;
	}

	public double[] encode(String t)
	{
		t = t.toUpperCase();
		int N = t.length();
		double ret[] = new double[3*N];
		int idx = 0;
		for(int i=0;i<t.length();i++)
		{
			char c = t.charAt(i);
			double z1 = (getZ1(c)-0.001923076923076976)/2.6160275521955336;
			double z2 = (getZ2(c)-0.0011538461538461635)/1.8589595518420015;
			double z3 = (getZ3(c)-0.0015384615384615096)/1.545268112160973;
			ret[idx]=z1; idx++;
			ret[idx]=z2; idx++;
			ret[idx]=z3; idx++;
		}
		return ret;
	}

	public static double getZ1(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return Z1[ ((char) x)-65 ];
	}	

	public static double getZ2(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return Z2[ ((char) x)-65 ];
	}

	public static double getZ3(char x)
	{
		if(x=='-'||x=='X')
			return 0.0;
		return Z3[ ((char) x)-65 ];
	}

}
