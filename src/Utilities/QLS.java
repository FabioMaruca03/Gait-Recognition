package Utilities;

public class QLS
{
	private static double a=(2+Math.sqrt(3));
	private static double b=(7-4*Math.sqrt(3));
	private static double xmax=20000;
	
	public static double F(int x)
	{
		double expB = Math.pow(b, (double)x/xmax);
		return (1-expB)/(a*expB + 1);
	}
}
