package Utilities;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.io.PrintWriter;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.Scanner;

import DataCenter.PersonWalk;
import DataCenter.Step;
//import ResultEvaluation.Result;
//import cwt.CWT;
//import cwt.CWT.Wavelet;
import flanagan.interpolation.CubicInterpolation;
import flanagan.interpolation.CubicSpline;
import flanagan.interpolation.LinearInterpolation;

public class Algorithm
{
	/**
	 * Basic formulation of DTW algorithm. NOTE: Works with int[].
	 * @return the DTW distance (as int value)
	 */
	public static int dynamicTimeWarping(int[] app1, int[] app2)
	{
		int DTW[][]=new int[app1.length+1][app2.length+1];
		//inizialization
		
		for(int i=1;i<DTW.length;i++)
			DTW[i][0]=Integer.MAX_VALUE;
		
		for(int i=1;i<DTW[0].length;i++)
			DTW[0][i]=Integer.MAX_VALUE;
		DTW[0][0]=0;
		
		//filling the matrix
		for(int i=1;i<DTW.length;i++)
		{
			for(int j=1;j<DTW[0].length;j++)
			{
				DTW[i][j]=(Math.abs(app1[i-1]-app2[j-1]));
				DTW[i][j]+=Integer.min(DTW[i-1][j], Integer.min(DTW[i][j-1],DTW[i-1][j-1]));
				//System.out.print(DTW[i][j]+" ");
			}
			//System.out.println();
		}
		return DTW[DTW.length-1][DTW[DTW.length-1].length-1];
	}
	
	
	/**
	 * Basic formulation of DTW algorithm. NOTE: Works with double[].
	 * @return the DTW distance (as double value)
	 */
	public static double dynamicTimeWarping(double[] app1, double[] app2)
	{
		double DTW[][]=new double[app1.length+1][app2.length+1];
		//inizialization
		
		for(int i=1;i<DTW.length;i++)
			DTW[i][0]=Integer.MAX_VALUE;
		
		for(int i=1;i<DTW[0].length;i++)
			DTW[0][i]=Integer.MAX_VALUE;
		DTW[0][0]=0;
		
		//filling the matrix
		for(int i=1;i<DTW.length;i++)
		{
			for(int j=1;j<DTW[0].length;j++)
			{
				DTW[i][j]=(Math.abs(app1[i-1]-app2[j-1]));
				DTW[i][j]+=Double.min(DTW[i-1][j], Double.min(DTW[i][j-1],DTW[i-1][j-1]));
				//System.out.print(DTW[i][j]+" ");
			}
			//System.out.println();
		}
		return DTW[DTW.length-1][DTW[DTW.length-1].length-1];
	}
	
	/**
	 * Basic formulation of DTW algorithm. NOTE: Works with int[].
	 * @return the DTW distance (as int value)
	 */
	public static int dynamicTimeWarpingSkyline(int[] app1, int[] app2, int wLen)
	{
		int DTW[][]=new int[app1.length+1][app2.length+1];
		//inizialization
		
		for(int i=1;i<DTW.length;i++)
			DTW[i][0]=Integer.MAX_VALUE;
		
		for(int i=1;i<DTW[0].length;i++)
			DTW[0][i]=Integer.MAX_VALUE;
		DTW[0][0]=0;
		
		//filling the matrix
		for(int i=1;i<DTW.length;i++)
		{
			for(int j=1;j<DTW[0].length;j++)
			{
				int iter=0;
				for(int k=Integer.max(0, i-wLen);k<Integer.min(app1.length, i+wLen);k++)
					for(int h=Integer.max(0, j-wLen);h<Integer.min(app2.length, j+wLen);h++)
					{
						DTW[i][j]+=(Math.abs(app1[k]-app2[h]));
						iter++;
					}
				DTW[i][j]/=iter;
				DTW[i][j]+=Integer.min(DTW[i-1][j], Integer.min(DTW[i][j-1],DTW[i-1][j-1]));
				//System.out.print(DTW[i][j]+" ");
			}
			//System.out.println();
		}
		return DTW[DTW.length-1][DTW[DTW.length-1].length-1];
	}
	
	/**
	 * Refined version of DTW. Takes as input the two walks to compare (and the related starts and stops computed during the walk creation - see also PersonWalk).<br>
	 * It applies the Z-Normalization (a.k.a. Gaussian normalization) to the gait signals.<br>
	 * It also uses the best endpoints adjustment in order to better fit the signals. See also https://www.cs.unm.edu/~mueen/DTW.pdf<br>
	 * WARNING the <b>c</b> must be chosen properly. It has to be proportional to the length of the walk to compare (and in general not greater than the 10% of the walk). 
	 * @param pw A PersonWalk to compare
	 * @param probe A PersonWalk to compare
	 * @param c window size for the search of endpoints
	 * @return DTW distance (as double)
	 */
//	public static double dynamicTimeWarpingZNormMultipe(PersonWalk pw, PersonWalk probe, int c)
//	{
//		int length1=pw.getLastStepEndIndex()-pw.getFirstStepIndex()+1;
//		int length2=probe.getLastStepEndIndex()-probe.getFirstStepIndex()+1;
//		double[] appX1=new double[length1], appY1=new double[length1], appZ1=new double[length1],
//				appX2=new double[length2], appY2=new double[length2], appZ2=new double[length2];
//		double DTW[][]=new double[length1+1][length2+1];
//		//inizialization
//		for(int i=pw.getFirstStepIndex();i<=Integer.min(pw.getLastStepEndIndex(), pw.getX().size());i++)
//		{
//			appX1[i-pw.getFirstStepIndex()]=pw.getX().get(i);
//			appY1[i-pw.getFirstStepIndex()]=pw.getY().get(i);
//			appZ1[i-pw.getFirstStepIndex()]=pw.getZ().get(i);
//		}
//		for(int i=probe.getFirstStepIndex();i<=Integer.min(probe.getLastStepEndIndex(), probe.getX().size()-1);i++)
//		{
//			appX2[i-probe.getFirstStepIndex()]=probe.getX().get(i);
//			appY2[i-probe.getFirstStepIndex()]=probe.getY().get(i);
//			appZ2[i-probe.getFirstStepIndex()]=probe.getZ().get(i);
//		}
//		//normalization
//		appX1=normalize(appX1);
//		appX2=normalize(appX2);
//		appY1=normalize(appY1);
//		appY2=normalize(appY2);
//		appZ1=normalize(appZ1);
//		appZ2=normalize(appZ2);
//		
//		for(int i=1;i<DTW.length;i++)
//			DTW[i][0]=Double.MAX_VALUE;
//		
//		for(int i=1;i<DTW[0].length;i++)
//			DTW[0][i]=Double.MAX_VALUE;
//		DTW[0][0]=0;
//			
//		//filling the matrix
//		for(int i=1;i<DTW.length;i++)
//		{
//			for(int j=1;j<DTW[0].length;j++)
//			{
//				DTW[i][j]=(Math.abs(appX1[i-1]-appX2[j-1])*Result.getWeightX()+
//						Math.abs(appY1[i-1]-appY2[j-1])*Result.getWeightY()+
//						Math.abs(appZ1[i-1]-appZ2[j-1])*Result.getWeightZ());
//				if(!((i==1 && j<=c) || (j==1 && i<=c)))
//					DTW[i][j]+=Double.min(DTW[i-1][j], Double.min(DTW[i][j-1],DTW[i-1][j-1]));
//			}
//		}
//
//		//endpoint chooser
//		double bestDis=DTW[DTW.length-1][DTW[0].length-1];
//		for(int i=DTW.length-c;i<DTW.length-1;i++)
//			if(DTW[i][DTW[0].length-1]<bestDis)
//				bestDis=DTW[i][DTW[0].length-1];
//		for(int j=DTW[0].length-c;j<DTW[0].length-1;j++)
//			if(DTW[DTW.length-1][j]<bestDis)
//				bestDis=DTW[DTW.length-1][j];
//			
//		return bestDis;
//	}
	
	/**
	 * Refined version of DTW. Takes as input the two steps to compare.<br>
	 * It applies the Z-Normalization (a.k.a. Gaussian normalization) to the gait signals.<br>
	 * It also uses the best endpoints adjustment in order to better fit the signals. See also https://www.cs.unm.edu/~mueen/DTW.pdf<br>
	 * WARNING the <b>c</b> must be chosen properly. It has to be proportional to the length of the walk to compare (and in general not greater than the 10% of the step). 
	 * @param pw A Step to compare
	 * @param probe A Step to compare
	 * @param c window size for the search of endpoints
	 * @return DTW distance (as double)
	 */
//	public static double dynamicTimeWarpingZNormMultipeStep(Step s, Step ss, int c)
//	{
//		int length1=s.getX().size();
//		int length2=ss.getX().size();
//		double[] appX1=new double[length1], appY1=new double[length1], appZ1=new double[length1],
//				appX2=new double[length2], appY2=new double[length2], appZ2=new double[length2];
//		double DTW[][]=new double[length1+1][length2+1];
//		//inizialization
//		for(int i=0;i<s.getX().size();i++)
//		{
//			appX1[i]=s.getX().get(i);
//			appY1[i]=s.getY().get(i);
//			appZ1[i]=s.getZ().get(i);
//		}
//		for(int i=0;i<ss.getX().size();i++)
//		{
//			appX2[i]=ss.getX().get(i);
//			appY2[i]=ss.getY().get(i);
//			appZ2[i]=ss.getZ().get(i);
//		}
//		//normalization
//		appX1=normalize(appX1);
//		appX2=normalize(appX2);
//		appY1=normalize(appY1);
//		appY2=normalize(appY2);
//		appZ1=normalize(appZ1);
//		appZ2=normalize(appZ2);
//		
//		for(int i=1;i<DTW.length;i++)
//			DTW[i][0]=Double.MAX_VALUE;
//		
//		for(int i=1;i<DTW[0].length;i++)
//			DTW[0][i]=Double.MAX_VALUE;
//		DTW[0][0]=0;
//		
//		//filling the matrix
//		for(int i=1;i<DTW.length;i++)
//		{
//			for(int j=1;j<DTW[0].length;j++)
//			{
//				DTW[i][j]=(Math.abs(appX1[i-1]-appX2[j-1])*Result.getWeightX()+
//						Math.abs(appY1[i-1]-appY2[j-1])*Result.getWeightY()+
//						Math.abs(appZ1[i-1]-appZ2[j-1])*Result.getWeightZ());
//				if(!((i==1 && j<=c) || (j==1 && i<=c)))
//					DTW[i][j]+=Double.min(DTW[i-1][j], Double.min(DTW[i][j-1],DTW[i-1][j-1]));
//			}
//		}
//		//endpoint chooser
//		double bestDis=DTW[DTW.length-1][DTW[0].length-1];
//		for(int i=DTW.length-c;i<DTW.length-1;i++)
//			if(DTW[i][DTW[0].length-1]<bestDis)
//				bestDis=DTW[i][DTW[0].length-1];
//		for(int j=DTW[0].length-c;j<DTW[0].length-1;j++)
//			if(DTW[DTW.length-1][j]<bestDis)
//				bestDis=DTW[DTW.length-1][j];
//			
//		return bestDis;
//	}
	
	/**
	 * Compute the magnitude vector from a PersonWalk
	 * @param pw PersonWalk at hand
	 * @return the magnitude of the walk
	 */
	public static double[] computeMagnitude(PersonWalk pw)
	{
		double[] magnitude=new double[pw.getY().size()];
		
		for(int i=0;i<pw.getY().size();i++)
			magnitude[i]=Math.sqrt(Math.pow(pw.getX().get(i), 2) + Math.pow(pw.getY().get(i), 2) + Math.pow(pw.getZ().get(i), 2));
		
		return magnitude;
	}
	
	/**
	 * Applies the Z-Normalization (a.k.a. Gaussian normalization) to a double[]<br>
	 * <b>NOTE: It does not override the content of the original array</b>
	 * @param toNormalize array to be normalized
	 * @return the normalized version of the array in input.
	 */
	public static double[] normalize(double[] toNormalize)
	{
		double mean=0;
		double sd=0;
		//mean
		for(double d : toNormalize)
			mean+=d;
		mean/=toNormalize.length;
		for(double d : toNormalize)
			sd+=(d-mean)*(d-mean);
		sd=Math.sqrt(sd/toNormalize.length);
		//Z normalization
		for(int i=0;i<toNormalize.length;i++)
			toNormalize[i]=(toNormalize[i]-mean)/sd;
		return toNormalize;
	}
	
	/**
	 * It computes the distance of two histograms.
	 * @param app1 An histogram to compare
	 * @param app2 The other histogram to compare
	 * @return histogram distance
	 */
	public static double histogramDistance(double[] app1, double[] app2)
	{
		//inizialization
		double dist=0;
		for(int i=1;i<app1.length;i++)
			dist+=Math.abs(app1[i]-app2[i]);
		return dist;
	}
	
	/**
	 * It creates the histograms for each axis starting from a PersonWalk
	 * NOTE: it only uses the portion of walk identified by [getFirstStepIndex(), getLastStepEndIndex()].
	 * @param pw The PersonWalk at hand
	 * @param bins # of the desired bins
	 * @return a [3 x bins] matrix in which rows represents the array of the histogram for the axes (order as x, y, and z) and the columns the value for the bin.
	 */ 
	public static double[][] createNormalizedHistogramFromWalk(PersonWalk pw, int bins, int min, int max)
	{
		return new double[][]{	createNormalizedHistrogramFromAxis(bins, pw.getX(), pw.getFirstStepIndex(), pw.getLastStepEndIndex(),max,min),
								createNormalizedHistrogramFromAxis(bins, pw.getY(), pw.getFirstStepIndex(), pw.getLastStepEndIndex(),max,min),
								createNormalizedHistrogramFromAxis(bins, pw.getZ(), pw.getFirstStepIndex(), pw.getLastStepEndIndex(),max,min)};	
	}
	
	/**
	 * It creates the histogram of a single axes. 
	 * Values in axes ranges from [-2,+2]
	 * @param bins # of desired bins
	 * @param axes the values for the axes at hand
	 * @param start starting point (used to select only a specific region or when you use the getFirstStepIndex() method of PersonWalk)
	 * @param end ending point (used to select only a specific region or when you use the getLastStepIndex() method of PersonWalk)
	 * @return An array containing the percentage of values in each bin (it is normalized by the size of the array)
	 */
	public static double[] createNormalizedHistrogramFromAxis(int bins, ArrayList<Integer> axes, int start, int end, int max, int min)
	{
		double[] histogram=new double[bins];
		double intervalLength=(max-min)*1.0/bins;
		for(int i=start; i<=(end<axes.size()? end : end-1); i++)
		{
			if(axes.get(i)<=min)
				histogram[0]++;
			else if(axes.get(i)>=max)
				histogram[bins-1]++;
			else
			{
				for(int j=1;j<=bins;j++)
					if(axes.get(i)<j*intervalLength)
					{
						histogram[j-1]++;
						break;
					}
			}
		}
		//Normalization for length
		for(int i=0; i<histogram.length;i++)
			histogram[i]/=axes.size();
		return histogram;
	}
	
	/**
	 * Compute the histogram distance from two PersonWalk.<br>
	 * The final score is given by distanceX * weightX + distanceY * weightY + distanceZ * weightZ<br>
	 * <b>NOTE It is suggested to have weightX, weightY, weightZ summing to 1000</b>
	 * @param pw
	 * @param probe
	 * @param bins
	 * @param weightX used for the combination for the final distance
	 * @param weightY used for the combination for the final distance
	 * @param weightZ used for the combination for the final distance
	 * @return
	 */
	public static double computeHistogramDistance(PersonWalk pw, PersonWalk probe, int bins, double weightX, double weightY, double weightZ, int max, int min)
	{
		double[][] histPw=createNormalizedHistogramFromWalk(pw, bins, max, min), histProbe=createNormalizedHistogramFromWalk(probe, bins, max, min);
		double disX=histogramDistance(histPw[0], histProbe[0]);
		double disY=histogramDistance(histPw[1], histProbe[1]);
		double disZ=histogramDistance(histPw[2], histProbe[2]);
//		System.out.println(disX+" "+disY+" "+disZ);
		return (Math.abs(disX*weightX)+Math.abs(disY*weightY)+Math.abs(disZ*weightZ));
	}
	
	/**
	 * A version of DTW in which the warping path is computed on the Y axis and then the SAME path is used in the comparison of the other two axes.
	 * It is worth pointing out that this does not guarantee the best comparison for the other two axis.
	 * @param pw gallery walk
	 * @param probe probe walk
	 * @param start if true use the firstStep and lastStep indexes (for WALK, otherwise the 0 and the size - SSW)
	 * @return distance values for X, Y, Z
	 */
	public static int[] dynamicTimeWarpingReconstruct(PersonWalk pw, PersonWalk probe, boolean start)
	{
		int[] app1=new int[pw.getY().size()], app2=new int[probe.getY().size()];
		if(start)
		{
			app1=new int[pw.getLastStepEndIndex()-pw.getFirstStepIndex()+1];
			app2=new int[probe.getLastStepEndIndex()-probe.getFirstStepIndex()+1];
			for(int i=pw.getFirstStepIndex(); i<(pw.getLastStepEndIndex()<pw.getY().size() ? pw.getLastStepEndIndex()+1 : pw.getY().size());i++)
				app1[i-pw.getFirstStepIndex()]=pw.getY().get(i);
			for(int i=probe.getFirstStepIndex(); i<(probe.getLastStepEndIndex()<probe.getY().size() ? probe.getLastStepEndIndex()+1 : probe.getY().size());i++)
				app2[i-probe.getFirstStepIndex()]=probe.getY().get(i);
		}
		else
		{
			for(int i=0; i<pw.getY().size();i++)
				app1[i]=pw.getY().get(i);
			for(int i=0; i<probe.getY().size();i++)
				app2[i]=probe.getY().get(i);
		}
		
		int DTW[][]=new int[app1.length+1][app2.length+1];
		//inizialization
		
		for(int i=1;i<DTW.length;i++)
			DTW[i][0]=Integer.MAX_VALUE;
		
		for(int i=1;i<DTW[0].length;i++)
			DTW[0][i]=Integer.MAX_VALUE;
		DTW[0][0]=0;
		
		//filling the matrix
		for(int i=1;i<DTW.length;i++)
		{
			for(int j=1;j<DTW[0].length;j++)
			{
				DTW[i][j]=(Math.abs(app1[i-1]-app2[j-1]));
				DTW[i][j]+=Integer.min(DTW[i-1][j], Integer.min(DTW[i][j-1],DTW[i-1][j-1]));
//				System.out.print(DTW[i][j]+"\t");
			}
//			System.out.println();
		}
		int distY=DTW[DTW.length-1][DTW[DTW.length-1].length-1];
		
		//reconstruct route
		ArrayList<Integer[]> route=new ArrayList<Integer[]>();
		int i=DTW.length-1, j=DTW[i].length-1;
		route.add(new Integer[]{i, j});
		while(i!=1 && j!=1)
		{
			Integer[] coord=new Integer[2];
			if(i>1 && j>1)
			{
				if(DTW[i-1][j-1]<=DTW[i][j-1] && DTW[i-1][j-1]<=DTW[i-1][j])
				{
					coord[0]=i-1;
					coord[1]=j-1;
					route.add(coord);
					i=i-1;
					j=j-1;
				}
				else if(DTW[i-1][j]<=DTW[i][j-1] && DTW[i-1][j]<=DTW[i-1][j-1])
				{
					coord[0]=i-1;
					coord[1]=j;
					route.add(coord);
					i=i-1;
				}
				else //if(DTW[i][j-1]<DTW[i-1][j] && DTW[i][j-1]<DTW[i-1][j-1])
				{
					coord[0]=i;
					coord[1]=j-1;
					route.add(coord);
					j=j-1;
				}
			}
			else if(i>1 && j==1)
			{
				coord[0]=i-1;
				coord[1]=j;
				route.add(coord);
				i=i-1;
			}
			else //if(i==i && j>1)
			{
				coord[0]=i;
				coord[1]=j-1;
				route.add(coord);
				j=j-1;
			}
		}
		
		int distX=0, distZ=0;
		for(Integer[] k : route)
		{
			distX+=Math.abs(pw.getX().get(k[0]-1)-probe.getX().get(k[1]-1));
			distZ+=Math.abs(pw.getZ().get(k[0]-1)-probe.getZ().get(k[1]-1));
		}
		
		return new int[]{distX, distY, distZ};
	}
	
	/*
	public static int cwtDistance(double[] app1, double[] app2)
	{
		ArrayList<Object> res1=computeCWT(app1);
		ArrayList<Object> res2=computeCWT(app2);
		
		return 0;
	}
	*/
	
	/*
	private static ArrayList<Object> computeCWT(double[] signal)
	{
		double dt=1;
		Wavelet mother=Wavelet.DOG;
		double param=2;
		double dj=0.25;
		double s0=0.00025;
		int jtot=49;
		return CWT.cWT(signal, dt, mother, param, s0, dj, jtot);
	}
	*/
	
	public static <T extends Number> ArrayList<Double> interpolatedNoiseRemoval(ArrayList<T> time, ArrayList<T> values, double dt, int future)
	{
		if(values==null || values.size()==0)
			return null;
		
		double threshold=0;
		for(T val : values)
			threshold+=val.doubleValue();
		threshold/=values.size();
		
		ArrayList<Double> interpolated=linearInterpolation(time, values, dt, null);
		interpolated=smoother((ArrayList)interpolated, (int)(0.05/dt), threshold);

//		ArrayList<Double>  interpolated=smoother((ArrayList)values, 25, 25);
//		interpolated=linearInterpolation(time, values, dt, null);
		
//		ArrayList<Double> interpolated=plateauInterpolation(time, values, dt, null);
//		ArrayList<Double> interpolated=plateauSmoothInterpolation(time, values, dt, null);
		ArrayList<Double> interDenoised=new ArrayList<Double>();
		
		if(values.size()>0)
			interDenoised.add(interpolated.get(0));
		//i=index of arraylist of values, k is the index of future
//		for(int i=1; i<interpolated.size();i++)
//		{
//			double current=interpolated.get(i).doubleValue();
//			double past=interpolated.get(i-1).doubleValue();
//			boolean desc=(current>past ? true : false);
//			boolean burst=false;
//			double avg=0;
//			int count=0;
//			for(int k=i+1; k*dt<=(Double)time.get(time.size()-1) && k<=i+1+future; k++)
//			{
//				double next=interpolated.get(k);
//				if(desc)
//				{
//					if(next>current)
//						burst=true;
//				}
//				else
//				{
//					if(next<current)
//						burst=true;
//				}
//				current=next;
//				avg+=next;
//				count++;
//			}
//			if(!burst)
//			{
////				interDenoised.add(interpolated.get(i));
////				System.out.println("!BURST");
//				for(int k=i+1; k*dt<=(Double)time.get(time.size()-1) && k<=i+1+(future/2); k++)
//					interDenoised.add(interpolated.get(k));
//				i+=(future/2);
//			}
//			else
//				interDenoised.add(count!=0 ? avg/count : 0);
//		}
		
		for(int i=1; i<interpolated.size();)
		{
			double current=interpolated.get(i).doubleValue();
			double past=interpolated.get(i-1).doubleValue();
			boolean desc=(current>past ? true : false);
			boolean burst=false;
			
			int k=i+1;
			for(; k<interpolated.size() && k<=i+1+future && !burst; k++)
			{
				double next=interpolated.get(k);
				if(desc)
				{
					if(next>current)
						burst=true;
				}
				else
				{
					if(next<current)
						burst=true;
				}
				current=next;
			}
			if(burst)
			{
				boolean newTrend=!desc;
				boolean newBurst=false;
				if(k<interpolated.size())
				{
					double pastFromNow=interpolated.get(k);
					int h=k+1;
					int burstIndex=h;
					for(; h<interpolated.size() && h<=k+1+future; h++)
					{
						double actFromNow=interpolated.get(h);
						if(newTrend)
						{
							if(actFromNow<pastFromNow)
							{
								newBurst=true;
								break;
							}
						}
						else
						{
							if(actFromNow>pastFromNow)
							{
								newBurst=true;
								break;
							}
						}
						pastFromNow=actFromNow;
					}
					if(newBurst)
					{
						if(h<interpolated.size())
						{
							double yA=interpolated.get(i);
							double yB=interpolated.get(h);
							double xA=dt*i;
							double xB=dt*h;
							
							double coeffAng=(yB-yA)/(xB-xA);
							double q=(xA*yB-xB*yA)/(xA-xB);
							
							interDenoised.add(interpolated.get(i));
							for(i=i+1; i<=h; i++)
								interDenoised.add(coeffAng*(dt*i)+q);
						}
					}
					else
					{
						if(i<burstIndex)
							for(; i<interpolated.size() && i<=burstIndex; i++)
								interDenoised.add(interpolated.get(i));
						else
						{
							interDenoised.add(interpolated.get(i));
							i++;
						}
					}
				}
				else
				{
					interDenoised.add(interpolated.get(i));
					i++;
				}
			}
			else
			{
				int start=i;
				for(; i<interpolated.size() && i<=start+(future/2); i++)
					interDenoised.add(interpolated.get(i));
			}
		}
		
		return interDenoised;
//		return interpolated;
	}
	
	public static ArrayList smoother(ArrayList<Number> y, int smoothLevel, double threshold)
	{
		if (y == null || y.size() == 0)
			return null;
		
		ArrayList<Double> smooth=new ArrayList<Double>();
				
		for(int i=0; i<smoothLevel; i++)
			smooth.add(y.get(i).doubleValue());
		
		for(int i=smoothLevel/2+1; i<y.size()-smoothLevel/2;i++)
		{
			double value=0;
			for(int j=i-smoothLevel/2;j<i+smoothLevel/2;j++)
				value+=y.get(j).doubleValue();
			value/=smoothLevel;
			if(Math.abs(value-y.get(i).doubleValue())<threshold)
				smooth.add(value);
			else
				smooth.add(y.get(i).doubleValue());
			
		}
		for(int i=y.size()-smoothLevel;i<y.size();i++)
			smooth.add(y.get(i).doubleValue());

		if (y.get(0) instanceof Integer)
		{
			ArrayList<Integer> retlist = new ArrayList<Integer>();
			for (Double d : smooth)
				retlist.add(d.intValue());
			return retlist;
		}
		return smooth;
	}
	
	public static ArrayList symmetricMovingWeightedAverage(ArrayList<Number> y, int N)
	{
		if (y == null || y.size() == 0)
			return null;
		
		ArrayList<Double> smooth=new ArrayList<Double>();
				
		for(int i=0; i<N/2+1; i++)
			smooth.add(y.get(i).doubleValue());
		
		for(int i=N/2+1; i<y.size()-N/2;i++)
		{
			double value=0;
			for(int j=i-N/2;j<i+N/2;j++)
				value+=y.get(j).doubleValue();
			value/=N;
			smooth.add(value);
		}
		for(int i=y.size()-N/2;i<y.size();i++)
			smooth.add(y.get(i).doubleValue());

		if (y.get(0) instanceof Integer)
		{
			ArrayList<Integer> retlist = new ArrayList<Integer>();
			for (Double d : smooth)
				retlist.add(d.intValue());
			return retlist;
		}
		return smooth;
	}
	
	public static ArrayList emphasize(ArrayList<Number> y, int emphaLevel, double threshold)
	{
		if (y == null || y.size() == 0)
			return null;
		
		ArrayList<Double> empha=new ArrayList<Double>();
		int i=0;
		
		for(;i<emphaLevel; i++)
			empha.add(y.get(i).doubleValue());
		
		for(; i<y.size()-emphaLevel;i++)
		{
			double sumDex=0, sumSix=0;
			for(int j=i-emphaLevel,h=i+emphaLevel;j<i;j++,h--)
			{
				sumDex+=y.get(h).doubleValue();
				sumSix+=y.get(j).doubleValue();
			}
			
			if(Math.abs(sumDex-sumSix)<threshold)
				empha.add(y.get(i).doubleValue()*2);
			else
				empha.add(y.get(i).doubleValue());
		}
		
		for(;i<y.size();i++)
			empha.add(y.get(i).doubleValue());
		
		if(y.get(0) instanceof Integer)
		{
			ArrayList<Integer> retlist = new ArrayList<Integer>();
			for (Double d : empha)
				retlist.add(d.intValue());
			return retlist;
		}
		
		return empha;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @param outFile if outFile is null, no file is created
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number> ArrayList<Double> linearInterpolation(ArrayList<T> x, ArrayList<T> y, double dt, String outFile)
	{
		ArrayList<Double> interpolated=new ArrayList<Double>();
		double appX[]=new double[x.size()+1];
		double appY[]=new double[y.size()+1];
		//the first element has to be in instant_time=0;
		appX[0]=0;
		appY[0]=y.get(0).doubleValue();
		for(int i=1; i<appX.length;i++)
		{
			appX[i]=x.get(i-1).doubleValue()-x.get(0).doubleValue(); //remove the start istant_time
			appY[i]=y.get(i-1).doubleValue();
		}
		
		//avoid prints of LinearInterpolation constructor
		PrintStream originalStream = System.out;
		PrintStream dummyStream = new PrintStream(new OutputStream()
		{
			public void write(int b) throws IOException{}
		});
		System.setOut(dummyStream);
//		
		LinearInterpolation.averageIdenticalAbscissae();
		LinearInterpolation li=new LinearInterpolation(appX, appY);
		
		System.setOut(originalStream);
		
		for(double i=0; i<=appX[appX.length-1]; i+=dt)
		{
			interpolated.add(li.interpolate(i));
		}
		
		if(outFile!=null)
			try
			{
				PrintWriter pw=new PrintWriter(outFile);
				for(int i=0; i<interpolated.size();i++)
					pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolated.get(i)+"\n");
				pw.close();
			} catch (FileNotFoundException e) {e.printStackTrace();}
		
		return interpolated;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number, V extends Number> ArrayList<Double> cubicInterpolation(ArrayList<T> x, ArrayList<T> y, double dt, int opt)
	{
		ArrayList<Double> interpolatedY=new ArrayList<Double>();
		double appX[]=new double[x.size()+1];
		double appY[]=new double[y.size()+1];
		appX[0]=0;
		appY[0]=y.get(0).doubleValue();
		for(int i=1; i<appX.length;i++)
		{
			appX[i]=x.get(i-1).doubleValue();
			appY[i]=y.get(i-1).doubleValue();
		}
		CubicInterpolation li=new CubicInterpolation(appX, appY, opt);
		
		for(double i=0; i<=appX[appX.length-1]; i+=dt)
		{
			interpolatedY.add(li.interpolate(i));
		}
		System.out.println(interpolatedY.size());
		try
		{
			PrintWriter pw=new PrintWriter("interpolatedCI"+opt+".csv");
			for(int i=0; i<interpolatedY.size();i++)
				pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolatedY.get(i)+"\n");
			pw.close();
		} catch (FileNotFoundException e) {e.printStackTrace();}
		
		return interpolatedY;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number> ArrayList<Double> splineInterpolation(ArrayList<T> x, ArrayList<T> y, double dt )
	{
		ArrayList<Double> interpolatedY=new ArrayList<Double>();
		double appX[]=new double[x.size()+1];
		double appY[]=new double[y.size()+1];
		appX[0]=0;
		appY[0]=y.get(0).doubleValue();
		for(int i=1; i<appX.length;i++)
		{
			appX[i]=x.get(i-1).doubleValue();
			appY[i]=y.get(i-1).doubleValue();
		}
		CubicSpline li=new CubicSpline(appX, appY);
		
		double[] app=appY.clone();
		Arrays.sort(app);
		double minValue=app[0];
		double maxValue=app[app.length-1];
		
		int j=0;
		for(double i=0; i<=appX[appX.length-1]; i+=dt, j++)
		{
			double interp=li.interpolate(i);
			if(interp<=maxValue && interp>=minValue)
				interpolatedY.add(interp);
			//if some round error occurs
			else
			{
				if(i>0)
					interpolatedY.add(interpolatedY.get(j-1));
				else
					interpolatedY.add(interpolatedY.get(0));
			}	
		}
		
		try
		{
			PrintWriter pw=new PrintWriter("interpolatedCS.csv");
			for(int i=0; i<interpolatedY.size();i++)
				pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolatedY.get(i)+"\n");
			pw.close();
		} catch (FileNotFoundException e) {e.printStackTrace();}
		
		return interpolatedY;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @param outFile if outFile is null, no file is created
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number> ArrayList<Double> plateauInterpolation(ArrayList<T> x, ArrayList<T> y, double dt, File outFile)
	{
		ArrayList<Double> interpolated=new ArrayList<Double>();
		if(y.size()>0)
		{
			interpolated.add((Double) y.get(0));
			int j=1, i=1;
			while(j*dt<=(Double)x.get(x.size()-1))
			{
				if(dt*j<(Double) x.get(i))
					interpolated.add((Double) y.get(i-1));
				else
				{
					interpolated.add((Double) y.get(i));
					i++;
				}
				j++;
			}
			
			if(outFile!=null)
				try
				{
					PrintWriter pw=new PrintWriter("interpolatedP.csv");
					for(i=0; i<interpolated.size();i++)
						pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolated.get(i)+"\n");
					pw.close();
				} catch (FileNotFoundException e) {e.printStackTrace();}
		}
		return interpolated;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @param outFile if outFile is null, no file is created
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number> ArrayList<Double> plateauSmoothInterpolation(ArrayList<T> x, ArrayList<T> y, double dt, File outFile)
	{
		ArrayList<Double> interpolated=new ArrayList<Double>();
		if(y.size()>0)
		{
			interpolated.add((Double) y.get(0));
			int j=1, i=1;
			while(j*dt<=(Double)x.get(x.size()-1))
			{
				if(Math.abs(dt*j-(Double) x.get(i-1))<Math.abs((Double) x.get(i)-dt*j))
					interpolated.add((Double) y.get(i-1));
				else
					interpolated.add((Double) y.get(i));
				if(dt*j>=(Double) x.get(i))
					i++;
				j++;
			}
			
			if(outFile!=null)
				try
				{
					PrintWriter pw=new PrintWriter("interpolatedPS.csv");
					for(i=0; i<interpolated.size();i++)
						pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolated.get(i)+"\n");
					pw.close();
				} catch (FileNotFoundException e) {e.printStackTrace();}
		}
		return interpolated;
	}
	
	/**
	 * 
	 * @param x time array
	 * @param y values array
	 * @param dt time distance desired
	 * @param outFile if outFile is null, no file is created
	 * @return the values array (each value is at <b>dt</b> distance)
	 */
	public static <T extends Number> ArrayList<Double> plateauSmoothP1DInterpolation(ArrayList<T> x, ArrayList<T> y, double dt, ArrayList<Integer> P1DMin, ArrayList<Integer> P1DMax, File outFile)
	{		
		ArrayList<Double> interpolated=new ArrayList<Double>();
		if(y.size()>0)
		{
			interpolated.add((Double) y.get(0));
			int j=1, i=1;
			while(j*dt<=(Double)x.get(x.size()-1))
			{
				if(P1DMin.contains(i) || P1DMax.contains(i))
					interpolated.add((Double) y.get(i));
				else
					if	(P1DMin.contains(i-1) || P1DMax.contains(i-1))
					interpolated.add((Double) y.get(i-1));
					//else the nearest
					else
						if(Math.abs(dt*j-(Double) x.get(i-1))<Math.abs((Double) x.get(i)-dt*j))
							interpolated.add((Double) y.get(i-1));
						else
							interpolated.add((Double) y.get(i));
				if(dt*j>=(Double) x.get(i))
					i++;
				j++;
			}
			
			if(outFile!=null)
				try
				{
					PrintWriter pw=new PrintWriter("interpolatedPS_P1D.csv");
					for(i=0; i<interpolated.size();i++)
						pw.write(String.format("%.3f",(dt*i)).replace(",", ".")+";"+interpolated.get(i)+"\n");
					pw.close();
				} catch (FileNotFoundException e) {e.printStackTrace();}
		}
		return interpolated;
	}
	
	/**
	 * 
	 * @param vec1 time*values vector1
	 * @param vec2 time*values vector2
	 * @return distance
	 */
	public static double dynamicTimeWarping2D(double[][] vec1, double[][] vec2)
	{
		double[][] normVec1=normalization(vec1);
		double[][] normVec2=normalization(vec2);
		
		double[][] distanceMatrix=new double[normVec1[0].length][normVec2[0].length];
		for(int i=0; i<normVec1[0].length;i++)
			for(int j=0; j<normVec2[0].length;j++)
//				distanceMatrix[i][j] = Math.abs(normVec1[0][i]-normVec2[0][j])+Math.abs(normVec1[1][i]-normVec2[1][j]); //base constraints
//				distanceMatrix[i][j] = Math.abs(normVec1[1][i]-normVec2[1][j]/Math.abs(normVec1[0][i]-normVec2[0][j]));
				distanceMatrix[i][j] = Math.abs(normVec1[1][i]-normVec2[1][j]*Math.abs(normVec1[0][i]-normVec2[0][j]));
		
		try
		{
			PrintWriter writer=new PrintWriter("distanceMatrix.csv");
			for(int i=0; i<distanceMatrix.length;i++)
			{
				for(int j=0; j<distanceMatrix[i].length;j++)
					writer.write(String.format("%.3f", distanceMatrix[i][j])+";");
				writer.write("\n");
			}
			writer.close();
		} catch (FileNotFoundException e) {e.printStackTrace();}
		
		//DTW Procedure
		double DTW[][]=new double[distanceMatrix.length+1][distanceMatrix[0].length+1];
		
		//inizialization
		for(int i=1;i<DTW.length;i++)
			DTW[i][0]=Double.MAX_VALUE;
		
		for(int i=1;i<DTW[0].length;i++)
			DTW[0][i]=Double.MAX_VALUE;
		DTW[0][0]=0;
		
		//filling the matrix
		for(int i=1;i<DTW.length;i++)
		{
			for(int j=1;j<DTW[0].length;j++)
			{
				DTW[i][j]=distanceMatrix[i-1][j-1];
				DTW[i][j]+=Double.min(DTW[i-1][j], Double.min(DTW[i][j-1], DTW[i-1][j-1]));
			}
		}
		
		try
		{
			PrintWriter writer=new PrintWriter("DTWMatrix.csv");
			for(int i=0; i<DTW.length;i++)
			{
				for(int j=0; j<DTW[i].length;j++)
					writer.write(String.format("%.3f", DTW[i][j])+";");
				writer.write("\n");
			}
			writer.close();
		} catch (FileNotFoundException e) {e.printStackTrace();}
		
		return DTW[DTW.length-1][DTW[DTW.length-1].length-1];
	}
	
	public static double[][] normalization(double[][] vec)
	{
		double meanTime=0, meanValues=0;
		
		for(int i=0; i<vec[0].length;i++)
		{
			meanTime+=vec[0][i];
			meanValues+=vec[1][i];
		}
		meanTime/=vec[0].length;
		meanValues/=vec[1].length;
		
		double stdTime=0, stdValues=0;
		for(int i=0; i<vec[0].length;i++)
		{
			stdTime+=((vec[0][i]-meanTime)*(vec[0][i]-meanTime));
			stdValues+=((vec[1][i]-meanValues)*(vec[1][i]-meanValues));
		}
		stdTime/=vec[0].length;
		stdValues/=vec[1].length;
		
		stdTime=Math.sqrt(stdTime);
		stdValues=Math.sqrt(stdValues);
		
		double[][] normalizedVec=new double[vec.length][vec[0].length];
		for(int i=0; i<normalizedVec[0].length;i++)
		{
			normalizedVec[0][i]=(vec[0][i]-meanTime)/stdTime;
			normalizedVec[1][i]=(vec[1][i]-meanValues)/stdValues;
		}
		return normalizedVec;
	}
	
//	public static void main(String[] args)
//	{
//		ArrayList<Double> sign=new ArrayList<Double>();
//		
//		try
//		{
//			Scanner in=new Scanner(new File("lucaM2.csv"));
//			while(in.hasNextLine())
//				sign.add(Double.parseDouble(in.nextLine().split(",")[1]));				
//			in.close();
//		}
//		catch(Exception e)
//		{
//			e.printStackTrace();
//		}
//		double signal[]=new double[sign.size()];
//		for(int i=0; i<sign.size();i++)
//			signal[i]=sign.get(i);
//		ArrayList<Object> app=computeCWT(signal);
//		
//		ComplexNumber[][] MRA=(ComplexNumber[][]) app.get(0);
//		double[] periods=(double[]) app.get(2);
//		double[] scales=(double[]) app.get(1);
//		double fourierFactor=(double) app.get(5);
//		System.out.println(fourierFactor);
//		
//		double[][] mraModulus=MatrixOps.modulus(MRA);
//		mraModulus=MatrixOps.transpose(mraModulus);
//		
//		double[] freq=new double[periods.length];
//		for(int i=0; i<freq.length; i++)
//		{
//			freq[i]=1/periods[i];
////			System.out.println(freq[i]);
//		}
//		
//			
//		XYSeriesCollection dataset=new XYSeriesCollection();
//		for(int i=0; i<mraModulus.length;i++)
//		{
//			XYSeries y=new XYSeries("Y"+i);
//			
//			for(int j=0; j<mraModulus[i].length;j++)
//			{
//				y.add(periods[i], mraModulus[i][j]);
////				y.add(freq[i], mraModulus[i][j]);
//			}
//			
//			dataset.addSeries(y);
//		}
//		
////		for(int i=0; i<MRA[0].length;i++)
////		{
////			XYSeries y=new XYSeries("Y"+i);
////			
////			for(int j=0; j<MRA[j].length;j++)
////			{
//////				y.add(freq[j], MRA[i][j].Real);
////				y.add(scales[j], MRA[j][i].Real);
////			}
////			
////			dataset.addSeries(y);
////		}
//		
//		JFreeChart y=ChartFactory.createXYLineChart("","", "", dataset, PlotOrientation.VERTICAL, false, true, false);
//		JFrame j=new JFrame("");
//		ChartPanel c=new ChartPanel(y);
//		j.add(c);
//		j.setVisible(true);
//		j.pack();
//		j.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		
//	}
	
	public static void zeroNorm(String dataset, String datasetOut, int mode)
	{
		new File(datasetOut).mkdir();
		for(File user : new File(dataset).listFiles())
		{
			File u=new File(datasetOut+File.separatorChar+user.getName());
			u.mkdir();
			for(File walk : user.listFiles())
			{
				File w=new File(datasetOut+File.separatorChar+user.getName()+File.separatorChar+walk.getName());
				w.mkdir();
				for(File walkN : walk.listFiles())
				{
					try
					{
						Scanner in = new Scanner(walkN);
					
						ArrayList<Double> time=new ArrayList<Double>();
						ArrayList<Double> appX=new ArrayList<Double>();
						ArrayList<Double> appY=new ArrayList<Double>();
						ArrayList<Double> appZ=new ArrayList<Double>();
						in.nextLine();
						while(in.hasNextLine())
						{
							String[] app=in.nextLine().replace(",",".").split(";",-1);
							time.add(Double.parseDouble(app[0]));
							appX.add(Double.parseDouble(app[1]));
							appY.add(Double.parseDouble(app[2]));
							appZ.add(Double.parseDouble(app[3]));
						}
						//average and sd
						double meanX=0, meanY=0, meanZ=0;
						double sdOverX=0, sdUnderX=0, sdOverY=0, sdUnderY=0, sdOverZ=0, sdUnderZ=0, sdX=0, sdY=0, sdZ=0;
						int countOverX=0, countUnderX=0, countOverY=0, countUnderY=0, countOverZ=0, countUnderZ=0;
						for(int i=0; i<appX.size();i++)
						{
							meanX+=appX.get(i);
							meanY+=appY.get(i);
							meanZ+=appZ.get(i);
						}
						meanX/=appX.size();
						meanY/=appY.size();
						meanZ/=appZ.size();
						for(int i=0; i<appX.size();i++)
						{
							if(appX.get(i)>=meanX)
							{
								sdOverX+=(appX.get(i)-meanX)*(appX.get(i)-meanX);
								countOverX++;
							}
							if(appX.get(i)<=meanX)
							{
								sdUnderX+=(appX.get(i)-meanX)*(appX.get(i)-meanX);
								countUnderX++;
							}
							if(appY.get(i)>=meanY)
							{
								sdOverY+=(appY.get(i)-meanY)*(appY.get(i)-meanY);
								countOverY++;
							}
							if(appY.get(i)<=meanY)
							{
								sdUnderY+=(appY.get(i)-meanY)*(appY.get(i)-meanY);
								countUnderY++;
							}
							if(appZ.get(i)>=meanZ)
							{
								sdOverZ+=(appZ.get(i)-meanZ)*(appZ.get(i)-meanZ);
								countOverZ++;
							}
							if(appZ.get(i)<=meanZ)
							{
								sdUnderZ+=(appZ.get(i)-meanZ)*(appZ.get(i)-meanZ);
								countUnderZ++;
							}
							sdX+=(appX.get(i)-meanX)*(appX.get(i)-meanX);
							sdY+=(appY.get(i)-meanY)*(appY.get(i)-meanY);
							sdZ+=(appZ.get(i)-meanZ)*(appZ.get(i)-meanZ);
						}
						sdOverX=Math.sqrt(sdOverX/countOverX);
						sdUnderX=Math.sqrt(sdUnderX/countUnderX);
						sdOverY=Math.sqrt(sdOverY/countOverY);
						sdUnderY=Math.sqrt(sdUnderY/countUnderY);
						sdOverZ=Math.sqrt(sdOverZ/countOverZ);
						sdUnderZ=Math.sqrt(sdUnderZ/countUnderZ);
						
						sdX=Math.sqrt(sdX/appX.size());
						sdY=Math.sqrt(sdY/appY.size());
						sdZ=Math.sqrt(sdZ/appZ.size());
						
						//Standardization
						ArrayList<Double> x=new ArrayList<Double>(), y=new ArrayList<Double>(), z=new ArrayList<Double>();
						for(int i=0;i<appX.size();i++)
						{
							switch(mode)
							{
								case 0://UpDownSD
									x.add(appX.get(i)>=meanX ? (appX.get(i)-meanX)/sdOverX : (appX.get(i)-meanX)/sdUnderX);
									y.add(appY.get(i)>=meanY ? (appY.get(i)-meanY)/sdOverY : (appY.get(i)-meanY)/sdUnderY);
									z.add(appZ.get(i)>=meanZ ? (appZ.get(i)-meanZ)/sdOverZ : (appZ.get(i)-meanZ)/sdUnderZ);
								case 1://NormSD
									x.add((appX.get(i)-meanX)/sdX);
									y.add((appY.get(i)-meanY)/sdY);
									z.add((appZ.get(i)-meanZ)/sdZ);
								case 2://WoSD
									x.add(appX.get(i)-meanX);
									y.add(appY.get(i)-meanY);
									z.add(appZ.get(i)-meanZ);
							}
						}
						PrintWriter pw=new PrintWriter(datasetOut+File.separatorChar+user.getName()+File.separatorChar+walk.getName()+File.separatorChar+walkN.getName());
						pw.write("time;x;y;z\n");
						for(int i=0; i<x.size();i++)
							pw.write(time.get(i)+";"+x.get(i)+";"+y.get(i)+";"+z.get(i)+"\n");
						pw.close();
					}
					catch (Exception e){e.printStackTrace();}
				}
			}
		}
	}
	
	public static void datasetSignalShaper(String dataset, String datasetOut)
	{
		for(File user : new File(dataset).listFiles())
		{
			for(File walk : user.listFiles())
			{
				for(File walkN : walk.listFiles())
				{
					try
					{
						Scanner in = new Scanner(walkN);
					
						ArrayList<Double> time=new ArrayList<Double>();
						ArrayList<Double> appX=new ArrayList<Double>();
						ArrayList<Double> appY=new ArrayList<Double>();
						ArrayList<Double> appZ=new ArrayList<Double>();
						in.nextLine();
						while(in.hasNextLine())
						{
							String[] app=in.nextLine().replace(",",".").split(";",-1);
							time.add(Double.parseDouble(app[0]));
							appX.add(Double.parseDouble(app[1]));
							appY.add(Double.parseDouble(app[2]));
							appZ.add(Double.parseDouble(app[3]));
						}
					}
					catch (Exception e){}
				}
			}
		}
	}
	
	public static ArrayList<Integer> signalShaper(ArrayList<Double> signal, ArrayList<Double> time)
	{
		ArrayList<Double> ret=new ArrayList<Double>();
		ArrayList<Double> ret2=new ArrayList<Double>();
		
//		int future=20;
//		if(signal.size()>0)
//		{
//			double last=signal.get(0);
//			
//			int start=0;
//			for(int i=1;i<signal.size();i++)
//			{
//				//if is a decreasing phase
//				if(last>=signal.get(i))
//				{
//					boolean decreasingPhase=true;
//					while(decreasingPhase && i<signal.size())
//					{
//						while(i<signal.size() && last>=signal.get(i))
//						{
//							last=signal.get(i);
//							i++;
//						}
//						if(i<signal.size())
//						{
//							boolean found=false;
//							for(int j=i;j<signal.size() && j<i+future && !found;j++)
//								if(last>signal.get(j))
//								{
//									last=signal.get(j);
//									i=j;
//									found=true;
//								}
//							if(!found)
//							{
//								System.out.println("inside "+start+" "+i);
//								i--;
//								//new signal
//								double yA=signal.get(start);
//								double yB=signal.get(i);
//								double xA=time.get(start);
//								double xB=time.get(i);
//								
//								double coeffAng=(yB-yA)/(xB-xA);
////								double q=(xA*yB-xB*yA)/(xA-xB);
//								double q=-coeffAng*xA+yA;
//								
//								for(int j=start; j<signal.size() && j<=i; j++)
//								{
////									System.out.print(String.format("%1$.3f\t", coeffAng*time.get(j)+q));
//									ret.add(coeffAng*time.get(j)+q);
//								}
//								
//								i++;
//								start=i;
//								last=signal.get(i);
//								decreasingPhase=false;
//							}
//						}
//					}
//				}
//				else
//				{
//					ret.add(signal.get(i));
//					start=i;
//				}
//				
//			}
//			ret2.add(ret.get(0));
//			for(int i=1;i<ret.size();i++)
//			{
//				int start=i;
//				//increasing phase
//				if(last>=ret.get(i))
//				{
//					last=ret.get(i);
//					while(i<ret.size() && last<=ret.get(i))
//					{
//						i++;
//						if(i<ret.size() && last>ret.get(i))
//						{
//							for(int j=i+1;j<ret.size() && j<i+future;j++)
//								if(ret.get(j)>last)
//								{
//									last=ret.get(j);
//									i=j+1;
//									break;
//								}
//						}
//						else
//						{
//							if(i<ret.size())
//								last=ret.get(i);
//							else
//								break;
//						}
//					}
//					//new signal
//					double yA=ret.get(start);
//					double yB=i<ret.size()? ret.get(i) : ret.get(ret.size()-1);
//					double xA=time.get(start);
//					double xB=i<time.size()? time.get(i) : time.get(time.size()-1);
//					
//					double coeffAng=(yB-yA)/(xB-xA);
//					double q=(xA*yB-xB*yA)/(xA-xB);
//					
//					ret2.add(ret.get(start));
//					for(int j=start+1; j<signal.size() && j<=i; j++)
//						ret2.add(coeffAng*time.get(j)+q);
//				}
//				else
//				{
//					ret2.add(ret.get(i));
//				}
//			}
//		}
		
		//time version
		double future=0.15;
		ArrayList<Integer> minIndexes=new ArrayList<Integer>();
		ArrayList<Integer> maxIndexes=new ArrayList<Integer>();
		//find all minimum
		if(signal.size()>0)
		{
			double last=signal.get(0);
			for(int i=1;i<signal.size();i++)
			{
				//if is a decreasing phase
				if(last>signal.get(i))
				{
					while(i<signal.size() && last>signal.get(i))
					{
						last=signal.get(i);
						i++;
					}
					if(i<signal.size())
					{
						minIndexes.add((i-1));
						last=signal.get(i);
					}
				}
				else
					last=signal.get(i);
			}
		}
		
		System.out.println("Min Ind");
		for(Integer k : minIndexes)
			System.out.print(k+"\t");
		System.out.println("size="+minIndexes.size()+"\n");
		
		//clear minimum of noise
		if(minIndexes.size()>0)
		{
			ArrayList<Integer> appMinimum=new ArrayList<Integer>();
			boolean change=true;
			while(change)
			{
				change=false;
				int last=minIndexes.get(0);
				for(int i=1; i<minIndexes.size();i++)
				{
					System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
					int act=minIndexes.get(i);
					if(Math.abs(time.get(last)-time.get(act))<future)
					{
						if(signal.get(last)<signal.get(act))
						{
							appMinimum.add(last);
							i++;
							change=true;
							if(i<minIndexes.size())
								last=minIndexes.get(i);
						}
						else
						{
							appMinimum.add(act);
							i++;
							if(i<minIndexes.size())
								last=minIndexes.get(i);
						}
					}
					else
					{
						appMinimum.add(last);
						last=minIndexes.get(i);
					}
				}
				
//				for(Integer k : minIndexes)
//					System.out.print(k+"\t");
//				System.out.println();
				for(Integer k : appMinimum)
					System.out.print(k+"\t");
				System.out.println("size="+appMinimum.size());
				minIndexes=new ArrayList<Integer>();
				minIndexes.addAll(appMinimum);
				appMinimum=new ArrayList<Integer>();
				
				//clear minimum of noise (from the end)
				if(minIndexes.size()>0)
				{
					last=minIndexes.get(minIndexes.size()-1);
					for(int i=minIndexes.size()-2; i>=0;i--)
					{
						System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
						int act=minIndexes.get(i);
						if(Math.abs(time.get(last)-time.get(act))<future)
						{
							if(signal.get(last)<signal.get(act))
							{
								appMinimum.add(last);
								i--;
								change=true;
								if(i>=0)
									last=minIndexes.get(i);
							}
							else
							{
								appMinimum.add(act);
								i--;
								if(i>=0)
									last=minIndexes.get(i);
							}
						}
						else
						{
							appMinimum.add(last);
							last=minIndexes.get(i);
						}
					}
						
//						for(Integer k : minIndexes)
//							System.out.print(k+"\t");
					System.out.println();
					for(Integer k : appMinimum)
						System.out.print(k+"\t");
					System.out.println("size="+appMinimum.size());
					minIndexes=new ArrayList<Integer>();
					minIndexes.addAll(appMinimum);
					appMinimum=new ArrayList<Integer>();
				}
			}
		}
		
//		return ret;
//		return ret2;
		Collections.sort(minIndexes);
		for(Integer k : minIndexes)
			System.out.println("index="+k+" time="+time.get(k));
		
		//find all maximum
		if(signal.size()>0)
		{
			double last=signal.get(0);
			for(int i=1;i<signal.size();i++)
			{
				//if is a decreasing phase
				if(last<signal.get(i))
				{
					while(i>signal.size() && last<signal.get(i))
					{
						last=signal.get(i);
						i++;
					}
					if(i<signal.size())
					{
						maxIndexes.add((i));
						last=signal.get(i);
					}
				}
				else
					last=signal.get(i);
			}
		}
		
		System.out.println("Max Ind");
		for(Integer k : maxIndexes)
			System.out.print(k+"\t");
		System.out.println("size="+maxIndexes.size()+"\n");
		
		//clear maximum of noise
		if(maxIndexes.size()>0)
		{
			ArrayList<Integer> appMaximum=new ArrayList<Integer>();
			boolean change=true;
			while(change)
			{
				change=false;
				int last=maxIndexes.get(0);
				for(int i=1; i<maxIndexes.size();i++)
				{
					System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
					int act=maxIndexes.get(i);
					if(Math.abs(time.get(last)-time.get(act))<future)
					{
						if(signal.get(last)>signal.get(act))
						{
							appMaximum.add(last);
							i++;
							change=true;
							if(i<maxIndexes.size())
								last=maxIndexes.get(i);
						}
						else
						{
							appMaximum.add(act);
							i++;
							if(i<maxIndexes.size())
								last=maxIndexes.get(i);
						}
					}
					else
					{
						appMaximum.add(last);
						last=maxIndexes.get(i);
					}
				}
				
//						for(Integer k : minIndexes)
//							System.out.print(k+"\t");
//						System.out.println();
				for(Integer k : appMaximum)
					System.out.print(k+"\t");
				System.out.println("size="+appMaximum.size());
				maxIndexes=new ArrayList<Integer>();
				maxIndexes.addAll(appMaximum);
				appMaximum=new ArrayList<Integer>();
				
				//clear maximum of noise (from the end)
				if(maxIndexes.size()>0)
				{
					last=maxIndexes.get(maxIndexes.size()-1);
					for(int i=maxIndexes.size()-2; i>=0;i--)
					{
						System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
						int act=maxIndexes.get(i);
						if(Math.abs(time.get(last)-time.get(act))<future)
						{
							if(signal.get(last)>signal.get(act))
							{
								appMaximum.add(last);
								i--;
								change=true;
								if(i>=0)
									last=maxIndexes.get(i);
							}
							else
							{
								appMaximum.add(act);
								i--;
								if(i>=0)
									last=maxIndexes.get(i);
							}
						}
						else
						{
							appMaximum.add(last);
							last=maxIndexes.get(i);
						}
					}
					
//								for(Integer k : minIndexes)
//									System.out.print(k+"\t");
					for(Integer k : appMaximum)
						System.out.print(k+"\t");
					System.out.println("size="+appMaximum.size());
					maxIndexes=new ArrayList<Integer>();
					maxIndexes.addAll(appMaximum);
					appMaximum=new ArrayList<Integer>();
				}
			}
		}
		
		Collections.sort(maxIndexes);
		for(Integer k : maxIndexes)
			System.out.println("MaxIndex="+k+" time="+time.get(k));
		
		//merge phase
		ArrayList<Integer> merge=new ArrayList<Integer>();
//		if(minIndexes.size()>0 && maxIndexes.size()>0)
//		{
//			if(minIndexes.get(0)<maxIndexes.get(0))
//			{
////				merge.add(minIndexes.get(0));
////				merge.add(maxIndexes.get(0));
////				int iMin=1, iMax=1;
////				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
////				{
////					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
////						iMin++;
////					if(iMin<minIndexes.size())
////						if(!merge.contains(minIndexes.get(iMin)))
////							merge.add(minIndexes.get(iMin));
////					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
////						iMax++;
////					if(iMax<maxIndexes.size())
////						if(!merge.contains(maxIndexes.get(iMax)))
////							merge.add(maxIndexes.get(iMax));
////				}
//				merge.add(minIndexes.get(0));
//				merge.add(maxIndexes.get(0));
//				int iMin=1, iMax=1;
//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//				{
//					int minI=iMin;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//					{
//						if(minIndexes.get(minI)<minIndexes.get(iMin))
//							minI=iMin;
//						iMin++;
//					}
//					if(iMin<minIndexes.size())
//						if(!merge.contains(minIndexes.get(minI)))
//							merge.add(minIndexes.get(minI));
//					int maxI=iMax;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//					{
//						if(maxIndexes.get(maxI)>maxIndexes.get(iMax))
//							maxI=iMax;
//						iMax++;
//					}
//					if(iMax<maxIndexes.size())
//						if(!merge.contains(maxIndexes.get(maxI)))
//							merge.add(maxIndexes.get(maxI));
//				}
//			}
//			else
//			{
////				merge.add(maxIndexes.get(0));
////				merge.add(minIndexes.get(0));
////				int iMin=1, iMax=1;
////				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
////				{
////					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
////						iMax++;
////					if(iMax<maxIndexes.size())
////						if(!merge.contains(maxIndexes.get(iMax)))
////							merge.add(maxIndexes.get(iMax));
////					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
////						iMin++;
////					if(iMin<minIndexes.size())
////						if(!merge.contains(minIndexes.get(iMin)))
////							merge.add(minIndexes.get(iMin));
////				}
//				
//				merge.add(maxIndexes.get(0));
//				merge.add(minIndexes.get(0));
//				int iMin=1, iMax=1;
//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//				{
//					int maxI=iMax;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//					{
//						if(maxIndexes.get(maxI)>maxIndexes.get(iMax))
//							maxI=iMax;
//						iMax++;
//					}
//					if(iMax<maxIndexes.size())
//						if(!merge.contains(maxIndexes.get(maxI)))
//							merge.add(maxIndexes.get(maxI));
//					int minI=iMin;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//					{
//						if(minIndexes.get(minI)<minIndexes.get(iMin))
//							minI=iMin;
//						iMin++;
//					}
//					if(iMin<minIndexes.size())
//						if(!merge.contains(minIndexes.get(minI)))
//							merge.add(minIndexes.get(minI));
//				}
//			}
//		}
		
		merge.addAll(minIndexes);
		merge.addAll(maxIndexes);
		Collections.sort(merge);
		
		System.out.println("\nMergeSize="+merge.size());
		for(Integer k : merge)
			System.out.println("mergedIndex="+k+" time="+time.get(k));
		
		PrintWriter pw;
		try
		{
			pw = new PrintWriter("min.csv");
			for(int i=0;i<minIndexes.size();i++)
				pw.write(minIndexes.get(i)+"\n");
			pw.close();
			pw = new PrintWriter("max.csv");
			for(int i=0;i<maxIndexes.size();i++)
				pw.write(maxIndexes.get(i)+"\n");
			pw.close();
			pw = new PrintWriter("merge.csv");
			for(int i=0;i<merge.size();i++)
				pw.write(merge.get(i)+"\n");
			pw.close();
		}
		catch (FileNotFoundException e)
		{
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return merge;
	}
	
	public static ArrayList<Double> signalShaperV4(ArrayList<Double> signal, ArrayList<Double> time)
	{
		ArrayList<Double> ret=new ArrayList<Double>();
		ArrayList<Double> ret2=new ArrayList<Double>();
		
		if(signal.size()>0)
		{
			double mean=0;
			double sd=0;
			for(Double d : signal)
				mean+=d;
			mean/=signal.size();
			for(Double d : signal)
				sd+=(d-mean)*(d-mean);
			sd=Math.sqrt(sd/signal.size());
			
			for(Double d : signal)
				ret.add(d-mean);
			
//			int smoothLevel=20;
//			double threshold=sd;
//			for(int i=0; i<smoothLevel/2+1; i++)
//				ret2.add(ret.get(i).doubleValue());
//			
//			for(int i=smoothLevel/2+1; i<ret.size()-smoothLevel/2;i++)
//			{
//				double value=0;
//				for(int j=i-smoothLevel/2;j<i+smoothLevel/2;j++)
//					value+=ret.get(j);
//				value/=smoothLevel;
//				if(Math.abs(value-ret.get(i))<threshold)
//				{
//					ret2.add(value);
//					System.out.println("here");
//				}
//				else
//					ret2.add(ret.get(i));
//				
//			}
//			for(int i=ret.size()-smoothLevel/2;i<ret.size();i++)
//				ret2.add(ret.get(i).doubleValue());

//			timeMode
			double timeSmoothLevel=0.20;
			double threshold=sd;
			int i=0;
			while(time.get(i)<timeSmoothLevel)
			{
				ret2.add(ret.get(i));
				i++;
			}
			
			while(time.get(time.size()-1)-time.get(i)>timeSmoothLevel)
			{
				//find start
				int start=i;
				while(start>0 && time.get(i)-time.get(start)<timeSmoothLevel)
					start--;
				int end=i;
				while(time.get(end)-time.get(i)<timeSmoothLevel)
					end++;
				double value=0;
				for(int j=start;j<=end;j++)
					value+=ret.get(j);
				value/=(end-start);
				System.out.println(value+" "+threshold);
				if(Math.abs(value-ret.get(i))<threshold)
				{
					ret2.add(value);
					System.out.println("here");
				}
				else
					ret2.add(ret.get(i));
				i++;
			}
			while(i<signal.size())
			{
				ret2.add(ret.get(i));
				i++;
			}
			
		}
		return ret2;
	}
	
	public static ArrayList<Integer> signalShaperV3(ArrayList<Double> signal, ArrayList<Double> time)
	{
		ArrayList<Double> ret=new ArrayList<Double>();
		ArrayList<Double> ret2=new ArrayList<Double>();
		
		ArrayList<Integer> merge=new ArrayList<Integer>();
		
		if(signal.size()>0)
		{
			double mean=0;
			double sdOver=0, sdUnder=0;
			int countOver=0, countUnder=0;
			for(Double d : signal)
				mean+=d;
			mean/=signal.size();
			for(Double d : signal)
			{
				if(d>=mean)
				{
					sdOver+=(d-mean)*(d-mean);
					countOver++;
				}
				if(d<=mean)
				{
					sdUnder+=(d-mean)*(d-mean);
					countUnder++;
				}
			}
			sdOver=Math.sqrt(sdOver/countOver);
			sdUnder=Math.sqrt(sdUnder/countUnder);
			
			System.out.println(mean+" "+sdUnder+" "+sdOver);
			
			//time version
			double future=0.10;
			ArrayList<Integer> minIndexes=new ArrayList<Integer>();
			ArrayList<Integer> maxIndexes=new ArrayList<Integer>();
			boolean min=true;
			//find all minimum
			if(signal.size()>0)
			{
				double last=signal.get(0);
				for(int i=1;i<signal.size();i++)
				{
					//if is a decreasing phase
					if(min)
					{
						while(i<signal.size() && (last>signal.get(i) || signal.get(i-1)<=mean-sdUnder))
						{
							last=signal.get(i);
							i++;
						}
						if(i<signal.size())
						{
							minIndexes.add((i-1));
							last=signal.get(i);
							min=false;
						}
					}
					else//if increasing phase
					{
						while(i<signal.size() && (last<signal.get(i) || signal.get(i-1)>=mean-sdOver))
						{
							last=signal.get(i);
							i++;
						}
						if(i<signal.size())
						{
							maxIndexes.add((i-1));
							last=signal.get(i);
							min=true;
						}
					}
				}
			}
			
			ArrayList<Integer> app=new ArrayList<>();
			app.addAll(minIndexes);
			app.addAll(maxIndexes);
			System.out.println(minIndexes.size()+" "+maxIndexes.size());
			for(Integer i:app)
				System.out.print(i+"\t");
			Collections.sort(app);
			System.out.println();
			for(Integer i:app)
				System.out.print(i+"\t");
			}
			//clear minimum of noise
//			if(minIndexes.size()>0)
//			{
//				ArrayList<Integer> appMinimum=new ArrayList<Integer>();
//				boolean change=true;
//				while(change)
//				{
//					change=false;
//					int last=minIndexes.get(0);
//					for(int i=1; i<minIndexes.size();i++)
//					{
//						System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
//						int act=minIndexes.get(i);
//						if(Math.abs(time.get(last)-time.get(act))<future)
//						{
//							if(signal.get(last)<signal.get(act))
//							{
//								appMinimum.add(last);
//								i++;
//								change=true;
//								if(i<minIndexes.size())
//									last=minIndexes.get(i);
//							}
//							else
//							{
//								appMinimum.add(act);
//								i++;
//								if(i<minIndexes.size())
//									last=minIndexes.get(i);
//							}
//						}
//						else
//						{
//							appMinimum.add(last);
//							last=minIndexes.get(i);
//						}
//						if(i==signal.size()-1)
//							appMinimum.add(minIndexes.get(i));
//					}
//					
//	//				for(Integer k : minIndexes)
//	//					System.out.print(k+"\t");
//	//				System.out.println();
//					for(Integer k : appMinimum)
//						System.out.print(k+"\t");
//					System.out.println("size="+appMinimum.size());
//					minIndexes=new ArrayList<Integer>();
//					minIndexes.addAll(appMinimum);
//					appMinimum=new ArrayList<Integer>();
//					
//					//clear minimum of noise (from the end)
//					if(minIndexes.size()>0)
//					{
//						last=minIndexes.get(minIndexes.size()-1);
//						for(int i=minIndexes.size()-2; i>=0;i--)
//						{
//							System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
//							int act=minIndexes.get(i);
//							if(Math.abs(time.get(last)-time.get(act))<future)
//							{
//								if(signal.get(last)<signal.get(act))
//								{
//									appMinimum.add(last);
//									i--;
//									change=true;
//									if(i>=0)
//										last=minIndexes.get(i);
//								}
//								else
//								{
//									appMinimum.add(act);
//									i--;
//									if(i>=0)
//										last=minIndexes.get(i);
//								}
//							}
//							else
//							{
//								appMinimum.add(last);
//								last=minIndexes.get(i);
//							}
//							if(i==0)
//								appMinimum.add(minIndexes.get(i));
//						}
//							
//	//						for(Integer k : minIndexes)
//	//							System.out.print(k+"\t");
//						System.out.println();
//						for(Integer k : appMinimum)
//							System.out.print(k+"\t");
//						System.out.println("size="+appMinimum.size());
//						minIndexes=new ArrayList<Integer>();
//						minIndexes.addAll(appMinimum);
//						appMinimum=new ArrayList<Integer>();
//					}
//				}
//			}
//			
//	//		return ret;
//	//		return ret2;
//			Collections.sort(minIndexes);
//			for(Integer k : minIndexes)
//				System.out.println("index="+k+" time="+time.get(k));
//			
//			//find all maximum
//			if(signal.size()>0)
//			{
//				double last=signal.get(0);
//				for(int i=1;i<signal.size();i++)
//				{
//					//if is a increasing phase
//					if(last<signal.get(i))
//					{
//						while(i<signal.size() && last<signal.get(i))
//						{
//							last=signal.get(i);
//							i++;
//						}
//						if(i<signal.size())
//						{
//							if(signal.get(i)>=mean+sdOver)
//								maxIndexes.add((i-1));
//							last=signal.get(i);
//						}
//					}
//					else
//						last=signal.get(i);
//				}
//			}
//			
//			System.out.println("Max Ind");
//			for(Integer k : maxIndexes)
//				System.out.print(k+"\t");
//			System.out.println("size="+maxIndexes.size()+"\n");
//			
//			//clear maximum of noise
//			if(maxIndexes.size()>0)
//			{
//				ArrayList<Integer> appMaximum=new ArrayList<Integer>();
//				boolean change=true;
//				while(change)
//				{
//					change=false;
//					int last=maxIndexes.get(0);
//					for(int i=1; i<maxIndexes.size();i++)
//					{
//						System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
//						int act=maxIndexes.get(i);
//						if(Math.abs(time.get(last)-time.get(act))<future)
//						{
//							if(signal.get(last)>signal.get(act))
//							{
//								appMaximum.add(last);
//								i++;
//								change=true;
//								if(i<maxIndexes.size())
//									last=maxIndexes.get(i);
//							}
//							else
//							{
//								appMaximum.add(act);
//								i++;
//								if(i<maxIndexes.size())
//									last=maxIndexes.get(i);
//							}
//						}
//						else
//						{
//							appMaximum.add(last);
//							last=maxIndexes.get(i);
//						}
//						if(i==appMaximum.size()-1)
//							appMaximum.add(maxIndexes.get(i));
//					}
//					
//	//						for(Integer k : minIndexes)
//	//							System.out.print(k+"\t");
//	//						System.out.println();
//					for(Integer k : appMaximum)
//						System.out.print(k+"\t");
//					System.out.println("size="+appMaximum.size());
//					maxIndexes=new ArrayList<Integer>();
//					maxIndexes.addAll(appMaximum);
//					appMaximum=new ArrayList<Integer>();
//					
//					//clear maximum of noise (from the end)
//					if(maxIndexes.size()>0)
//					{
//						last=maxIndexes.get(maxIndexes.size()-1);
//						for(int i=maxIndexes.size()-2; i>=0;i--)
//						{
//							System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
//							int act=maxIndexes.get(i);
//							if(Math.abs(time.get(last)-time.get(act))<future)
//							{
//								if(signal.get(last)>signal.get(act))
//								{
//									appMaximum.add(last);
//									i--;
//									change=true;
//									if(i>=0)
//										last=maxIndexes.get(i);
//								}
//								else
//								{
//									appMaximum.add(act);
//									i--;
//									if(i>=0)
//										last=maxIndexes.get(i);
//								}
//							}
//							else
//							{
//								appMaximum.add(last);
//								last=maxIndexes.get(i);
//							}
//							if(i==0)
//								appMaximum.add(maxIndexes.get(i));
//						}
//						
//	//								for(Integer k : minIndexes)
//	//									System.out.print(k+"\t");
//						for(Integer k : appMaximum)
//							System.out.print(k+"\t");
//						System.out.println("size="+appMaximum.size());
//						maxIndexes=new ArrayList<Integer>();
//						maxIndexes.addAll(appMaximum);
//						appMaximum=new ArrayList<Integer>();
//					}
//				}
//			}
//			
//			Collections.sort(maxIndexes);
//			for(Integer k : maxIndexes)
//				System.out.println("MaxIndex="+k+" time="+time.get(k));
//			
//			//merge phase
//			if(minIndexes.size()>0 && maxIndexes.size()>0)
//			{
//				if(minIndexes.get(0)<maxIndexes.get(0))
//				{
//	//				merge.add(minIndexes.get(0));
//	//				merge.add(maxIndexes.get(0));
//	//				int iMin=1, iMax=1;
//	//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//	//				{
//	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//	//						iMin++;
//	//					if(iMin<minIndexes.size())
//	//						if(!merge.contains(minIndexes.get(iMin)))
//	//							merge.add(minIndexes.get(iMin));
//	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//	//						iMax++;
//	//					if(iMax<maxIndexes.size())
//	//						if(!merge.contains(maxIndexes.get(iMax)))
//	//							merge.add(maxIndexes.get(iMax));
//	//				}
//					merge.add(minIndexes.get(0));
//					merge.add(maxIndexes.get(0));
//					int iMin=1, iMax=1;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//					{
//						int minI=iMin;
//						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//						{
//							if(minIndexes.get(iMin)<minIndexes.get(minI))
//								minI=iMin;
//							iMin++;
//						}
//						if(!merge.contains(minIndexes.get(minI)))
//								merge.add(iMin<=minIndexes.size() ? minIndexes.get(minI) : minIndexes.get(minIndexes.size()-1));
//						
//						int maxI=iMax;
//						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//						{
//							if(maxIndexes.get(iMax)>maxIndexes.get(maxI))
//								maxI=iMax;
//							iMax++;
//						}
//						if(!merge.contains(maxIndexes.get(maxI)))
//								merge.add(iMax<maxIndexes.size() ? maxIndexes.get(maxI) : maxIndexes.get(maxIndexes.size()-1));
//					}
//				}
//				else
//				{
//	//				merge.add(maxIndexes.get(0));
//	//				merge.add(minIndexes.get(0));
//	//				int iMin=1, iMax=1;
//	//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//	//				{
//	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//	//						iMax++;
//	//					if(iMax<maxIndexes.size())
//	//						if(!merge.contains(maxIndexes.get(iMax)))
//	//							merge.add(maxIndexes.get(iMax));
//	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//	//						iMin++;
//	//					if(iMin<minIndexes.size())
//	//						if(!merge.contains(minIndexes.get(iMin)))
//	//							merge.add(minIndexes.get(iMin));
//	//				}
//					
//					merge.add(maxIndexes.get(0));
//					merge.add(minIndexes.get(0));
//					int iMin=1, iMax=1;
//					while(iMin<minIndexes.size() && iMax<maxIndexes.size())
//					{
//						int maxI=iMax;
//						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
//						{
//							if(maxIndexes.get(maxI)>maxIndexes.get(iMax))
//								maxI=iMax;
//							iMax++;
//						}
//						if(iMax<maxIndexes.size())
//						{
//							if(!merge.contains(maxIndexes.get(maxI)))
//								merge.add(maxIndexes.get(maxI));
//						}
//						else
//						{
//							if(!merge.contains(maxIndexes.get(maxI)))
//								merge.add(maxIndexes.get(maxIndexes.size()-1));
//						}
//						
//						int minI=iMin;
//						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
//						{
//							if(minIndexes.get(iMin)<minIndexes.get(minI))
//								minI=iMin;
//							iMin++;
//						}
//						if(!merge.contains(minIndexes.get(minI)))
//								merge.add(iMin<=minIndexes.size() ? minIndexes.get(minI) : minIndexes.get(minIndexes.size()-1));
//					}
//				}
//			}
//			
////			merge.addAll(minIndexes);
////			merge.addAll(maxIndexes);
////			Collections.sort(merge);
//			
//			System.out.println("\nMergeSize="+merge.size());
//			for(Integer k : merge)
//				System.out.println("mergedIndex="+k+" time="+time.get(k));
//			
//			PrintWriter pw;
//			try
//			{
//				pw = new PrintWriter("min.csv");
//				for(int i=0;i<minIndexes.size();i++)
//					pw.write(minIndexes.get(i)+"\n");
//				pw.close();
//				pw = new PrintWriter("max.csv");
//				for(int i=0;i<maxIndexes.size();i++)
//					pw.write(maxIndexes.get(i)+"\n");
//				pw.close();
//				pw = new PrintWriter("merge.csv");
//				for(int i=0;i<merge.size();i++)
//					pw.write(merge.get(i)+"\n");
//				pw.close();
//			}
//			catch (FileNotFoundException e)
//			{
//				// TODO Auto-generated catch block
//				e.printStackTrace();
//			}
//		}
//		return merge;
		return null;
	}
	
	public static ArrayList<Integer> signalShaperV2(ArrayList<Double> signal, ArrayList<Double> time)
	{
		ArrayList<Double> ret=new ArrayList<Double>();
		ArrayList<Double> ret2=new ArrayList<Double>();
		
		ArrayList<Integer> merge=new ArrayList<Integer>();
		
		if(signal.size()>0)
		{
			double mean=0;
			double sdOver=0, sdUnder=0;
			int countOver=0, countUnder=0;
			for(Double d : signal)
				mean+=d;
			mean/=signal.size();
			for(Double d : signal)
			{
				if(d>=mean)
				{
					sdOver+=(d-mean)*(d-mean);
					countOver++;
				}
				if(d<=mean)
				{
					sdUnder+=(d-mean)*(d-mean);
					countUnder++;
				}
			}
			sdOver=Math.sqrt(sdOver/countOver);
			sdUnder=Math.sqrt(sdUnder/countUnder);
			
			System.out.println(mean+" "+sdUnder+" "+sdOver);
			
			//time version
			double future=0.10;
			ArrayList<Integer> minIndexes=new ArrayList<Integer>();
			ArrayList<Integer> maxIndexes=new ArrayList<Integer>();
			//find all minimum
			if(signal.size()>0)
			{
				double last=signal.get(0);
				for(int i=1;i<signal.size();i++)
				{
					//if is a decreasing phase
					if(last>signal.get(i))
					{
						while(i<signal.size() && last>signal.get(i))
						{
							last=signal.get(i);
							i++;
						}
						if(i<signal.size())
						{
							if(signal.get(i-1)<=mean-sdUnder)
								minIndexes.add((i-1));
							last=signal.get(i);
						}
					}
					else
						last=signal.get(i);
				}
			}
			
			System.out.println("Min Ind");
			for(Integer k : minIndexes)
				System.out.print(k+"\t");
			System.out.println("size="+minIndexes.size()+"\n");
			
			//clear minimum of noise
			if(minIndexes.size()>0)
			{
				ArrayList<Integer> appMinimum=new ArrayList<Integer>();
				boolean change=true;
				while(change)
				{
					change=false;
					int last=minIndexes.get(0);
					for(int i=1; i<minIndexes.size();i++)
					{
						System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
						int act=minIndexes.get(i);
						if(Math.abs(time.get(last)-time.get(act))<future)
						{
							if(signal.get(last)<signal.get(act))
							{
								appMinimum.add(last);
								i++;
								change=true;
								if(i<minIndexes.size())
									last=minIndexes.get(i);
							}
							else
							{
								appMinimum.add(act);
								i++;
								if(i<minIndexes.size())
									last=minIndexes.get(i);
							}
						}
						else
						{
							appMinimum.add(last);
							last=minIndexes.get(i);
						}
						if(i==signal.size()-1)
							appMinimum.add(minIndexes.get(i));
					}
					
	//				for(Integer k : minIndexes)
	//					System.out.print(k+"\t");
	//				System.out.println();
					for(Integer k : appMinimum)
						System.out.print(k+"\t");
					System.out.println("size="+appMinimum.size());
					minIndexes=new ArrayList<Integer>();
					minIndexes.addAll(appMinimum);
					appMinimum=new ArrayList<Integer>();
					
					//clear minimum of noise (from the end)
					if(minIndexes.size()>0)
					{
						last=minIndexes.get(minIndexes.size()-1);
						for(int i=minIndexes.size()-2; i>=0;i--)
						{
							System.out.println("last="+last +" minI="+ minIndexes.get(i)+" i="+i);
							int act=minIndexes.get(i);
							if(Math.abs(time.get(last)-time.get(act))<future)
							{
								if(signal.get(last)<signal.get(act))
								{
									appMinimum.add(last);
									i--;
									change=true;
									if(i>=0)
										last=minIndexes.get(i);
								}
								else
								{
									appMinimum.add(act);
									i--;
									if(i>=0)
										last=minIndexes.get(i);
								}
							}
							else
							{
								appMinimum.add(last);
								last=minIndexes.get(i);
							}
							if(i==0)
								appMinimum.add(minIndexes.get(i));
						}
							
	//						for(Integer k : minIndexes)
	//							System.out.print(k+"\t");
						System.out.println();
						for(Integer k : appMinimum)
							System.out.print(k+"\t");
						System.out.println("size="+appMinimum.size());
						minIndexes=new ArrayList<Integer>();
						minIndexes.addAll(appMinimum);
						appMinimum=new ArrayList<Integer>();
					}
				}
			}
			
	//		return ret;
	//		return ret2;
			Collections.sort(minIndexes);
			for(Integer k : minIndexes)
				System.out.println("index="+k+" time="+time.get(k));
			
			//find all maximum
			if(signal.size()>0)
			{
				double last=signal.get(0);
				for(int i=1;i<signal.size();i++)
				{
					//if is a increasing phase
					if(last<signal.get(i))
					{
						while(i<signal.size() && last<signal.get(i))
						{
							last=signal.get(i);
							i++;
						}
						if(i<signal.size())
						{
							if(signal.get(i)>=mean+sdOver)
								maxIndexes.add((i-1));
							last=signal.get(i);
						}
					}
					else
						last=signal.get(i);
				}
			}
			
			System.out.println("Max Ind");
			for(Integer k : maxIndexes)
				System.out.print(k+"\t");
			System.out.println("size="+maxIndexes.size()+"\n");
			
			//clear maximum of noise
			if(maxIndexes.size()>0)
			{
				ArrayList<Integer> appMaximum=new ArrayList<Integer>();
				boolean change=true;
				while(change)
				{
					change=false;
					int last=maxIndexes.get(0);
					for(int i=1; i<maxIndexes.size();i++)
					{
						System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
						int act=maxIndexes.get(i);
						if(Math.abs(time.get(last)-time.get(act))<future)
						{
							if(signal.get(last)>signal.get(act))
							{
								appMaximum.add(last);
								i++;
								change=true;
								if(i<maxIndexes.size())
									last=maxIndexes.get(i);
							}
							else
							{
								appMaximum.add(act);
								i++;
								if(i<maxIndexes.size())
									last=maxIndexes.get(i);
							}
						}
						else
						{
							appMaximum.add(last);
							last=maxIndexes.get(i);
						}
						if(i==appMaximum.size()-1)
							appMaximum.add(maxIndexes.get(i));
					}
					
	//						for(Integer k : minIndexes)
	//							System.out.print(k+"\t");
	//						System.out.println();
					for(Integer k : appMaximum)
						System.out.print(k+"\t");
					System.out.println("size="+appMaximum.size());
					maxIndexes=new ArrayList<Integer>();
					maxIndexes.addAll(appMaximum);
					appMaximum=new ArrayList<Integer>();
					
					//clear maximum of noise (from the end)
					if(maxIndexes.size()>0)
					{
						last=maxIndexes.get(maxIndexes.size()-1);
						for(int i=maxIndexes.size()-2; i>=0;i--)
						{
							System.out.println("last="+last +" maxI="+ maxIndexes.get(i)+" i="+i);
							int act=maxIndexes.get(i);
							if(Math.abs(time.get(last)-time.get(act))<future)
							{
								if(signal.get(last)>signal.get(act))
								{
									appMaximum.add(last);
									i--;
									change=true;
									if(i>=0)
										last=maxIndexes.get(i);
								}
								else
								{
									appMaximum.add(act);
									i--;
									if(i>=0)
										last=maxIndexes.get(i);
								}
							}
							else
							{
								appMaximum.add(last);
								last=maxIndexes.get(i);
							}
							if(i==0)
								appMaximum.add(maxIndexes.get(i));
						}
						
	//								for(Integer k : minIndexes)
	//									System.out.print(k+"\t");
						for(Integer k : appMaximum)
							System.out.print(k+"\t");
						System.out.println("size="+appMaximum.size());
						maxIndexes=new ArrayList<Integer>();
						maxIndexes.addAll(appMaximum);
						appMaximum=new ArrayList<Integer>();
					}
				}
			}
			
			Collections.sort(maxIndexes);
			for(Integer k : maxIndexes)
				System.out.println("MaxIndex="+k+" time="+time.get(k));
			
			//merge phase
			if(minIndexes.size()>0 && maxIndexes.size()>0)
			{
				if(minIndexes.get(0)<maxIndexes.get(0))
				{
	//				merge.add(minIndexes.get(0));
	//				merge.add(maxIndexes.get(0));
	//				int iMin=1, iMax=1;
	//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
	//				{
	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
	//						iMin++;
	//					if(iMin<minIndexes.size())
	//						if(!merge.contains(minIndexes.get(iMin)))
	//							merge.add(minIndexes.get(iMin));
	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
	//						iMax++;
	//					if(iMax<maxIndexes.size())
	//						if(!merge.contains(maxIndexes.get(iMax)))
	//							merge.add(maxIndexes.get(iMax));
	//				}
					merge.add(minIndexes.get(0));
					merge.add(maxIndexes.get(0));
					int iMin=1, iMax=1;
					while(iMin<minIndexes.size() && iMax<maxIndexes.size())
					{
						int minI=iMin;
						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
						{
							if(minIndexes.get(iMin)<minIndexes.get(minI))
								minI=iMin;
							iMin++;
						}
						if(!merge.contains(minIndexes.get(minI)))
								merge.add(iMin<=minIndexes.size() ? minIndexes.get(minI) : minIndexes.get(minIndexes.size()-1));
						
						int maxI=iMax;
						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
						{
							if(maxIndexes.get(iMax)>maxIndexes.get(maxI))
								maxI=iMax;
							iMax++;
						}
						if(!merge.contains(maxIndexes.get(maxI)))
								merge.add(iMax<maxIndexes.size() ? maxIndexes.get(maxI) : maxIndexes.get(maxIndexes.size()-1));
					}
				}
				else
				{
	//				merge.add(maxIndexes.get(0));
	//				merge.add(minIndexes.get(0));
	//				int iMin=1, iMax=1;
	//				while(iMin<minIndexes.size() && iMax<maxIndexes.size())
	//				{
	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
	//						iMax++;
	//					if(iMax<maxIndexes.size())
	//						if(!merge.contains(maxIndexes.get(iMax)))
	//							merge.add(maxIndexes.get(iMax));
	//					while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
	//						iMin++;
	//					if(iMin<minIndexes.size())
	//						if(!merge.contains(minIndexes.get(iMin)))
	//							merge.add(minIndexes.get(iMin));
	//				}
					
					merge.add(maxIndexes.get(0));
					merge.add(minIndexes.get(0));
					int iMin=1, iMax=1;
					while(iMin<minIndexes.size() && iMax<maxIndexes.size())
					{
						int maxI=iMax;
						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && maxIndexes.get(iMax)<minIndexes.get(iMin))
						{
							if(maxIndexes.get(maxI)>maxIndexes.get(iMax))
								maxI=iMax;
							iMax++;
						}
						if(iMax<maxIndexes.size())
						{
							if(!merge.contains(maxIndexes.get(maxI)))
								merge.add(maxIndexes.get(maxI));
						}
						else
						{
							if(!merge.contains(maxIndexes.get(maxI)))
								merge.add(maxIndexes.get(maxIndexes.size()-1));
						}
						
						int minI=iMin;
						while(iMin<minIndexes.size() && iMax<maxIndexes.size() && minIndexes.get(iMin)<maxIndexes.get(iMax))
						{
							if(minIndexes.get(iMin)<minIndexes.get(minI))
								minI=iMin;
							iMin++;
						}
						if(!merge.contains(minIndexes.get(minI)))
								merge.add(iMin<=minIndexes.size() ? minIndexes.get(minI) : minIndexes.get(minIndexes.size()-1));
					}
				}
			}
			
//			merge.addAll(minIndexes);
//			merge.addAll(maxIndexes);
//			Collections.sort(merge);
			
			System.out.println("\nMergeSize="+merge.size());
			for(Integer k : merge)
				System.out.println("mergedIndex="+k+" time="+time.get(k));
			
			PrintWriter pw;
			try
			{
				pw = new PrintWriter("min.csv");
				for(int i=0;i<minIndexes.size();i++)
					pw.write(minIndexes.get(i)+"\n");
				pw.close();
				pw = new PrintWriter("max.csv");
				for(int i=0;i<maxIndexes.size();i++)
					pw.write(maxIndexes.get(i)+"\n");
				pw.close();
				pw = new PrintWriter("merge.csv");
				for(int i=0;i<merge.size();i++)
					pw.write(merge.get(i)+"\n");
				pw.close();
			}
			catch (FileNotFoundException e)
			{
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
		}
		return merge;
	}

	
	public static void convertDatasetNameAndPositionsForFeatureExtraction(File dataIn, File dataOut)
	{
		dataOut.mkdirs();
		for(File f : dataIn.listFiles())
		{
			String subjID=f.getName();
			for(File ff : f.listFiles())
			{
				String nWalk=ff.getName().replaceAll("[A-Za-z]", "");
				for(File fff : ff.listFiles())
				{
					try
					{
						String content=new String(Files.readAllBytes(Paths.get(fff.getAbsolutePath())), StandardCharsets.UTF_8);
						PrintWriter writer=new PrintWriter(dataOut.getAbsolutePath()+File.separatorChar+subjID+"_"+nWalk+".csv");
						writer.write(content);
						writer.close();
					}
					catch (IOException e)
					{
						e.printStackTrace();
					}
				}
			}
		}
	}
	
	public static void selecterWalkNForFeatureExtraction(ArrayList<Integer> positions, File dataIn, File dataOut)
	{
		dataOut.mkdirs();
		for(File f : dataIn.listFiles())
		{
			String filename=f.getName();
			if(positions.contains(Integer.parseInt(filename.split("_")[2].replaceAll(".csv", ""))))
				try
				{
					String content=new String(Files.readAllBytes(Paths.get(f.getAbsolutePath())), StandardCharsets.UTF_8);
					PrintWriter writer=new PrintWriter(dataOut.getAbsolutePath()+File.separatorChar+filename);
					writer.write(content);
					writer.close();
				}
				catch (IOException e)
				{
					e.printStackTrace();
				}
		}
	}
	
	public static void main(String[] args) throws FileNotFoundException{}
	
}
