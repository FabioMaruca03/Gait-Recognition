package DataCenter;

import static AppClasses.Constants.*;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintWriter;
import java.lang.reflect.Array;
import java.sql.Timestamp;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Scanner;

import AppClasses.MaximumPos;
//import StepCounter.SignalShower;
//import StepCounter.StepCounter;
import Utilities.Algorithm;

public class PersonWalk
{
	
	private ArrayList<Double> arrayDoubleInterpolationX;
	private ArrayList<Double> arrayDoubleInterpolationY;
	private ArrayList<Double> arrayDoubleInterpolationZ;
	
	//WALK DATA
	//Signal data
	private ArrayList<Integer> x=new ArrayList<Integer>();
	private ArrayList<Integer> y=new ArrayList<Integer>();
	private ArrayList<Integer> z=new ArrayList<Integer>();
	private ArrayList<Integer> time=new ArrayList<Integer>();
	
	//Computed data - These data are setted to default values, but are then recomputed on demanding.
	private int threshold1step=105;			//threshold for the first step start. Setted to a default value
	private int thresholdLastStepEnd=110;	//threshold for the last step end. Setted to a default value
	private int firstStepIndex=0;					//index of the first step start
	private int lastStepEndIndex;					//index of the last step end
	private int stepsPerWalk=-1;				//default number of walk steps.
	private double mean;							//mean of values of y axis
	private double standardDeviation;				//standard deviation of values y axis
	/*
	 * The first time mean and standardDeviation are computed on the entire signal. These values are used to compute "threshold1step"
	 * and "thresholdLastStepEnd" and to detect the "firstStepIndex" and the "lastStepEndIndex". After that are recomputed on the signal
	 * part between the two indexes.
	*/
	//File details
	private File originalFile;						//original walk file reference
	//Device Info
	private String deviceInfo;						//Model and Brand of the mobile device
	private Timestamp timestamp;					//Timestamp del file
	//Walk Acquisition Information
	public enum states{TRUE,FALSE,NA};				//general enum for the data acquisition conditions. NA stands for unknown			
	private states beacon;							//beacon used, not used or unknown
	private states normalized;						//normalization used, not used or unknown
	
	//After Segmentation
	private ArrayList<Step> steps=new ArrayList<Step>();
	private ArrayList<Step> stepsWithOutlier=new ArrayList<Step>();
	
	private Step bestStep=null;
	
	//Increase segmentation accuracy variables
	
	int distanceAfter=20; //after a maxima, indicates the number of following signal samples to be analyzed before
						  //accept that maxima as a step start/stop point.
	
	//Segmentation Variables (to compute)
	private int stepThreshold; //stepThreshold value used in segmentation algorithm
	private int stepEquilibrium; // stepEquilibrium value used in segmentation algorithm
	
	//Other variables
	private static int smoothLevel=10; //interval length for smoothing function
	private static int smoothThreshold=25; //threshold for smoothing function
	
	private ArrayList<Integer> P1DMaxIndexes=new ArrayList<Integer>();
	private ArrayList<Integer> P1DMinIndexes=new ArrayList<Integer>();
	
	/*public PersonWalk(ArrayList<Integer> x, ArrayList<Integer> y, ArrayList<Integer> z, ArrayList<Integer> time)
	{
		this.x=x;
		this.y=y;
		this.z=z;
		this.time=time;
	}*/
	
	/**
	 * Create a PersonWalk without any data except the steps and entire walk.
	 * @param steps step composing the walk
	 */
	public PersonWalk(ArrayList<Step> steps)
	{
		firstStepIndex=0;
		for(Step s: steps)
		{
			time.addAll(s.getTime());
			x.addAll(s.getX());
			y.addAll(s.getY());
			z.addAll(s.getZ());
		}
		lastStepEndIndex=x.size()-1;
	}
	
	/**
	 * Create a PersonWalk based on <b>pw</b> signal.
	 * @param pw base signal
	 * @param startIndex start of new walk on <b>pw</b> signal
	 * @param endIndex end of new walk on <b>pw</b> signal
	 */
	public PersonWalk(PersonWalk pw, int startIndex, int endIndex)
	{
		firstStepIndex=startIndex;
		for(int i=startIndex; i<(endIndex<pw.getY().size() ? endIndex : pw.getY().size()-1);i++)
		{
			time.add(pw.getTime().get(i));
			x.add(pw.getX().get(i));
			y.add(pw.getY().get(i));
			z.add(pw.getZ().get(i));
		}
		lastStepEndIndex=endIndex;
		P1DMaxIndexes=pw.P1DMaxIndexes;
		P1DMinIndexes=pw.P1DMinIndexes;
	}
	
	/**
	 * Create a PersonWalk with data in the passed file
	 * @param walk file with signal/user information
	 */
	public PersonWalk(File walk)
	{
		if(walk.isFile())
		{
				
			try
			{
				Scanner in=new Scanner(walk);
				in.nextLine();
				while(in.hasNextLine())
				{
					String s[]=in.nextLine().split(";",-1);
					
					if(s!=null && s.length==4)
					{					
						s[0]=s[0].replace(",", ".");
						s[1]=s[1].replace(",", ".");
						s[2]=s[2].replace(",", ".");
						s[3]=s[3].replace(",", ".");
						
						Float timeApp=Float.parseFloat(s[0])*1000;
						Float xApp=Float.parseFloat(s[1])*100;
						Float yApp=Float.parseFloat(s[2])*100;
						Float zApp=Float.parseFloat(s[3])*100;
						
						time.add(new Integer(timeApp.intValue()));
						x.add(new Integer(xApp.intValue()));
						y.add(new Integer(yApp.intValue()));
						z.add(new Integer(zApp.intValue()));
					}
					else
					{
						if(s.length==1)
						{
							try
							{
								//the string will be in the format: "Total Steps:x"
								stepsPerWalk=Integer.parseInt(s[0].split(":",-1)[1]);
								System.out.println(stepsPerWalk);
								if(in.hasNextLine())
								{
									s=in.nextLine().split(":",-1);
									if(s.length==2)
									{
										s=s[1].split("_");
										deviceInfo=s[0];
										timestamp=new Timestamp(Integer.parseInt(s[2]), Integer.parseInt(s[3]), Integer.parseInt(s[4]), Integer.parseInt(s[5]), Integer.parseInt(s[6]), Integer.parseInt(s[7]), 0);
									}
								}
							}
							catch(Exception e)
							{
								stepsPerWalk=10;
							}
						}
					}
				}
				in.close();
			}
			catch (Exception e)
			{
//				e.printStackTrace();
				System.out.println(walk);
			}
			originalFile=walk;
			
//			File fP1Dmax=new File("temp"+File.separatorChar+"indexes"+File.separatorChar+getOriginalFile().getName().substring(0, getOriginalFile().getName().length()-4)+"YMaximaIndexes.txt");
//			File fP1Dmin=new File("temp"+File.separatorChar+"indexes"+File.separatorChar+getOriginalFile().getName().substring(0, getOriginalFile().getName().length()-4)+"YMinimaIndexes.txt");
//			
//			try
//			{
//				if(!fP1Dmax.exists())
//				{
//					File outname=new File("temp"+File.separatorChar+"onlyY"+File.separatorChar+getOriginalFile().getName().substring(0, getOriginalFile().getName().length()-4)+"Y.txt");
//					PrintWriter yAxisW=new PrintWriter(outname);
//					for(double i : getY())
//						yAxisW.write((i/100)+"\n");
//					yAxisW.close();
//					Runtime rt = Runtime.getRuntime();
//					Process proc = rt.exec("filterPoints.exe "+outname+" "+"temp"+File.separatorChar+"indexes"+File.separatorChar);
//					proc.waitFor();
//				}
//				
//				Scanner in=new Scanner(fP1Dmax);
//				while(in.hasNextLine())
//					P1DMaxIndexes.add(Integer.parseInt(in.nextLine()));
//				
//				if(!fP1Dmin.exists())
//				{
//					File outname=new File("temp"+File.separatorChar+"onlyY"+File.separatorChar+getOriginalFile().getName().substring(0, getOriginalFile().getName().length()-4)+"Y.txt");
//					if(!outname.exists())
//					{
//						PrintWriter yAxisW=new PrintWriter(outname);
//						for(double i : getY())
//							yAxisW.write((i/100)+"\n");
//						yAxisW.close();
//					}
//					
//					Runtime rt = Runtime.getRuntime();
//					Process proc = rt.exec("filterPoints.exe "+outname+" "+"temp"+File.separatorChar+"indexes"+File.separatorChar);
//					proc.waitFor();
//				}
//				
//				in=new Scanner(fP1Dmin);
//				while(in.hasNextLine())
//					P1DMinIndexes.add(Integer.parseInt(in.nextLine()));
//				in.close();
//			}catch (Exception e) {e.printStackTrace();}
			
//			double dt=2;
//			ArrayList<Double> xInt=Algorithm.linearInterpolation(time, x, dt, null);
//			ArrayList<Double> yInt=Algorithm.linearInterpolation(time, y, dt, null);
//			ArrayList<Double> zInt=Algorithm.linearInterpolation(time, z, dt, null);	
//			x=new ArrayList<Integer>();
//			y=new ArrayList<Integer>();
//			z=new ArrayList<Integer>();
//			time=new ArrayList<Integer>();
//			for(int i=0; i<xInt.size(); i++)
//			{
//				time.add(i*(int)dt);
//				x.add(xInt.get(i).intValue());
//				y.add(yInt.get(i).intValue());
//				z.add(zInt.get(i).intValue());
//			}
			
//			try
//			{
//				PrintWriter pw=new PrintWriter(walk+"int");
//				for(int i=0; i<y.size();i++)
//				{
//					pw.write(time.get(i)+";"+x.get(i)+";"+y.get(i)+";"+z.get(i)+"\n");
//				}
//				pw.close();
//			}catch (Exception e) {
//				// TODO: handle exception
//			}
			
			lastStepEndIndex=y.size()-1;
			
			computeMeanAndStandardDeviation();

			threshold1step=new Double(mean+(standardDeviation)).intValue();
			thresholdLastStepEnd=new Double(mean+(standardDeviation)).intValue();
			
			findAndSetFirstStep(); //TODO test with different values
			findAndSetLastStepEnd(); //TODO test with different values
			
			firstStepIndex=0;
			lastStepEndIndex=x.size()-1;
			
			computeMeanAndStandardDeviation();

			ArrayList<Integer> acc[]=new ArrayList[3];
//			acc[0]=x;
//			acc[1]=y;
//			acc[2]=z;
//			
			//if stepsPerWalk is still the default -1 means that the file does not contains the last two rows, so it must be computed.
//			if(stepsPerWalk==-1)
//				stepsPerWalk=StepCounter.stepCounter(acc);
			//System.out.println(stepsPerWalk);
			
//			stepsPerWalk=10;
			
			//collapse maxima
//			computeP1DModV2(P1DMaxIndexes);
//			stepsPerWalk=P1DMaxIndexes.size()-1;
//			System.out.println(stepsPerWalk);
			stepsPerWalk=10;
		}
	}
	
	/**
	 * Do nothing constructor
	 */
	public PersonWalk()
	{
		
	}

	/**
	 * It computes the mean and the standard deviation values for the y axis and then assign them to the relative variables
	 */
	public void computeMeanAndStandardDeviation()
	{
		mean=0;
		for(Integer i : y)
			mean+=i;
		mean/=y.size();
		for(Integer i : y)
			standardDeviation+=(i-mean)*(i-mean);
		standardDeviation/=y.size();
		standardDeviation=Math.sqrt(standardDeviation);
	}
	
	/**
	 * It finds the last point of useful part of signal. It is based on the <b>thresholdLastStep</b> variable. Used for initial and final part of signal noise reduction.
	 */
	public void findAndSetLastStepEnd()
	{
		if(y.size()>0)
		{
			boolean found=false;
			int i=y.size()-1;
			for(;i>=0 && !found;i--)
				if(y.get(i)>thresholdLastStepEnd)
					found=true;
			lastStepEndIndex=i;
			for(i=i-1;i>=0;i--)
			{
				if(y.get(i)>=y.get(i+1))
					lastStepEndIndex=i;
				else
					return;
			}
		}
	}
	
	/**
	 * It finds the first point of useful part of signal. It is based on the <b>threshold1step</b> variable. Used for initial and final part of signal noise reduction.
	 */
	public void findAndSetFirstStep()
	{
		if(y.size()>0)
		{
			boolean found=false;
			int i=0;
			for(;i<y.size() && !found;i++)
				if(y.get(i)>threshold1step)
					found=true;
			firstStepIndex=i;
			for(i=i+1;i<y.size();i++)
			{
				if(y.get(i)>=y.get(i-1))
					firstStepIndex=i;
				else
					return;
			}
		}
	}
	
	/*
	private void computeP1DModV2(ArrayList<Integer> maximaIndexes)
	{
		ArrayList<MaximumPos> goodMaxima=new ArrayList<MaximumPos>();
		int range=25*2;
		for(int i=1;i<maximaIndexes.size();i++)
		{
			int act=maximaIndexes.get(i);
			int prec=maximaIndexes.get(i-1);
			int next=(i==maximaIndexes.size()-1) ? -1 : maximaIndexes.get(i+1);
			
			if(act-prec>=range)
			{
				//2.1 - 2.4 - 1.1 - 1.4 - 3.1 - 3.4 - 4.1 - 4.4
				if(next!=-1)
					goodMaxima.add(new MaximumPos(y.get(prec), prec));	
				else
				{
					goodMaxima.add(new MaximumPos(y.get(prec), prec));
					goodMaxima.add(new MaximumPos(y.get(act), act));
				}
			}
			else
			{
				//act-prec<range
				if(y.get(prec)<=y.get(act))
				{
					//1.3 - 1.2 - 2.2 - 2.3
					if(next!=-1)
					{
						if(next-act<=range)
						{
							//1.3 - 2.3
							if(y.get(next)<=y.get(act))
							{
								//1.3
								goodMaxima.add(new MaximumPos(y.get(act), act));
								if(i+2<maximaIndexes.size())
									i++;
							}
							else
							{
								//2.3
								if(goodMaxima.size()>0)
								{
									MaximumPos last=goodMaxima.get(goodMaxima.size()-1);
									//if the last inserted is distant, add the present prec
									if(prec-last.getPos()>range)
									{
										goodMaxima.add(new MaximumPos(y.get(prec), prec));
									}
									else
										;//do nothing and go further
								}
								else
								;//do nothing and go further
							}
						}
						else
						{
							//next-act>50
							//1.2 - 2.2
							goodMaxima.add(new MaximumPos(y.get(act), act));
							if(i+2<maximaIndexes.size())
								i++;
						}
					}
					else
					{
						goodMaxima.add(new MaximumPos(y.get(act), act));
					}
				}
				else
				{
					//prec>act
					//4.2 - 4.3 - 3.2 - 3.3
					if(next!=-1)
					{
						if(y.get(next)<=y.get(act))
						{
							//4.2 - 4.3
							goodMaxima.add(new MaximumPos(y.get(prec), prec));
						}
						else
						{
							//3.3 - 3.2
							if(next-act<range)
							{
								//3.3
								if(goodMaxima.size()>0)
								{
									MaximumPos last=goodMaxima.get(goodMaxima.size()-1);
									if(prec-last.getPos()>range)
									{
										goodMaxima.add(new MaximumPos(y.get(prec), prec));
										if(next-prec<range)
										{
											i+=2;
										}
									}
									else
									{
										if(last.getValue()<y.get(prec))
										{
											goodMaxima.remove(last);
											goodMaxima.add(new MaximumPos(y.get(prec), prec));
										}
										//else do nothing and go further
									}
								}
							}
							else
							{
								//3.2
								if(goodMaxima.size()>0)
								{
									MaximumPos last=goodMaxima.get(goodMaxima.size()-1);
									if(prec-last.getPos()>range)
									{
										goodMaxima.add(new MaximumPos(y.get(prec), prec));
									}
									else
									{
										goodMaxima.remove(last);
										goodMaxima.add(new MaximumPos(y.get(prec), prec));
									}
								}
								else
								{
									goodMaxima.add(new MaximumPos(y.get(prec), prec));
								}
							}
						}
					}
					else
					{
						goodMaxima.add(new MaximumPos(y.get(prec), prec));
						break;
					}
				}
			}
		}
		maximaIndexes=new ArrayList<Integer>();
		for(MaximumPos mp : goodMaxima)
			if(!maximaIndexes.contains(mp.getPos()))
				maximaIndexes.add(mp.getPos());
	}
	*/
	
	/**
	 * String representation of a walk
	 */
	public String toString()
	{
		String walk=originalFile+"\n";
		walk+="x: ";
		for(Integer i:x)
			walk+=i+" ";
		walk+="\ny: ";
		for(Integer i:y)
			walk+=i+" ";
		walk+="\nz: ";
		for(Integer i:z)
			walk+=i+" ";
		walk+="\n";
		walk+="First Step Index = "+firstStepIndex+" First Step Value = "+y.get(firstStepIndex)+"\n";
		walk+="Last Step End Index = "+lastStepEndIndex+" Last Step End Value = "+y.get(lastStepEndIndex)+"\n";
		return walk;
	}
	
	/**
	 * Create the file with start and stop points of the segmented steps 
	 * @param f start folder. The file will be created in a sub-folder (by default called "data"). 
	 */
	public void createSteps(File f)
	{
		File dataDir=new File(f.getAbsolutePath()+"\\"+extrapolatedDataFolderName);
		
		if(!dataDir.exists()) {
			if(dataDir.mkdir())
			{
				separateAndWriteSteps(dataDir);
			}
			else
			{
				System.out.println("Error in DATA folder creation !!!!!!!");
				System.out.println(f.getAbsolutePath());
			}
		}
	}
	
	/**
	 * Execute all procedure of step segmentation algorithm and write all files in the related folder
	 * @param dataDir destination folder.
	 */
	private void separateAndWriteSteps(File dataDir)
	{
        stepEquilibrium=(int)(mean-standardDeviation);
        
		chooseStepThreshold();

		//separate steps
		separateSteps(this);

		//delete outlier
		deleteOutlierAndChooseTheBest(this);
		
		//write steps
		printSteps(dataDir);
		
		//create file with stepThreshold //create file with the best step
		try
		{
			PrintWriter p=new PrintWriter(new File(dataDir+"\\"+stepThresholdFileName));
			p.write(stepThreshold+"");
			p.close();
			
			p=new PrintWriter(new File(dataDir+"\\"+bestStepFileName));
			p.write(bestStep.toString());
			p.close();
			
			p=new PrintWriter(new File(dataDir+"\\"+stepEquilibriumFileName));
			p.write(stepEquilibrium+"");
			p.close();
		}
		catch (Exception e)
		{
			System.out.println(originalFile);
			e.printStackTrace();
		}
	}

	public void chooseStepThreshold()
	{
        ArrayList<MaximumPos> relativeMax=new ArrayList<MaximumPos>();

        for(int i=0;i<y.size();i++)
        {
        	int last=y.get(i);
            i++;
            
            //looks for a minima
            while(i<=lastStepEndIndex && i<y.size())
            {
            	if(y.get(i)>last)
                {
            		last=y.get(i);
                    break;
                }
                else
                {
                    last=y.get(i);
                    i++;
                }
            }
            
            //looks for a maxima over the stepEquilibrium value
            while(i<=lastStepEndIndex && i<y.size())
            {
                if(y.get(i)>stepEquilibrium && y.get(i)<last)
                {
                	//look if is a real possible interesting Maximum
                    int k=i+1;
                    MaximumPos max=new MaximumPos(y.get(i), i);

                    //looks if is the real interesting maxima
                    while(k<i+distanceAfter && k<y.size())
                    {
                        if(y.get(k)>max.getValue())
                            max=new MaximumPos(y.get(k), k);
                        k++;
                    }
                    
                    if(max.getPos()!=i)
                    {
                        relativeMax.add(max);
                        i=max.getPos();
                    }
                    else
                        relativeMax.add(new MaximumPos(y.get(i-1), i-1));
                    break;
                }
                else
                {
                	last=y.get(i);
                    i++;
                }
            }
        }
		Collections.sort(relativeMax);
		//assuming the PersonWalks are made of 10 steps each one
		int k=stepsPerWalk;
		stepThreshold=(relativeMax.size()>=k) ?	relativeMax.get(relativeMax.size()-k).getValue() : relativeMax.get(0).getValue();
//		System.out.println("Theshold "+stepThreshold);
	}
	
	/**
	 * ORIGINAL VERSION
	 * It performs the step segmentation of signal based on stepEquilibrium and stepThreshold variables 
	 * @param pw walk to be segmented
	 */
	public static void separateSteps(PersonWalk pw)
	{
		pw.steps=new ArrayList<Step>();
		pw.stepsWithOutlier=new ArrayList<Step>();
                
        //find first step start
        int i=pw.getFirstStepIndex();
		while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
		{
			if(pw.getY().get(i)>=pw.stepThreshold)
				break;
            else
            	i++;
		}
        int app=0;
        if(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
            app=pw.getY().get(i);
        else
            return;
        for(i=i+1;i<=pw.getLastStepEndIndex() && i<pw.getY().size();)
        {
            if(pw.getY().get(i)<app)
                break;
            else
            {
                app=pw.getY().get(i);
                i++;
            }
        }
        i--;
        //looks if it is the real interesting maxima
        MaximumPos max=new MaximumPos(pw.y.get(i), i);
        for(int k=i;k<i+pw.distanceAfter && k<=pw.lastStepEndIndex && k<pw.y.size();)
        {
            if(pw.y.get(k)>max.getValue())
                max=new MaximumPos(pw.y.get(k), k);
            k++;
        }
        if(max.getPos()!=i)
            i=max.getPos();
        
        //finds all other steps
		while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
		{
			int firstMaxIndex=i;
			
			//find stepEquilibrium after the step begin
			while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
			{
				if(pw.getY().get(i)<pw.stepEquilibrium)
					break;
				else
					i++;
			}
                        
			//find signal point over stepThreshold
			while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
			{
				if(pw.getY().get(i)>=pw.stepThreshold)
					break;
				else
					i++;
			}
			//find maxima
			if(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
            {
                app=pw.getY().get(i);
                for(i=i+1;i<=pw.getLastStepEndIndex() && i<pw.getY().size();)
                {
                    if(pw.getY().get(i)<app)
                    	break;
                    else
                    {
                        app=pw.getY().get(i);
                        i++;
                    }
                }
            }
			
			//looks if it is the real interesting maxima
			if(i<pw.y.size())
			{
		        max=new MaximumPos(pw.y.get(i), i);
		        for(int k=i;k<i+pw.distanceAfter && k<=pw.lastStepEndIndex && k<pw.y.size();)
		        {
		            if(pw.y.get(k)>max.getValue())
		                max=new MaximumPos(pw.y.get(k), k);
		            k++;
		        }
		        if(max.getPos()!=i)
		            i=max.getPos();
	            
		        //step creation
				if(i>pw.getLastStepEndIndex())
					break;
				else
	            	pw.steps.add(new Step(pw,firstMaxIndex, i-1));
				i--;
				
				//avoid loops
	            if(firstMaxIndex==i)
	            	break;
			}
			else
				break;
		}
	}
	
	//WITH OVERLAP
//	public static void separateSteps(PersonWalk pw)
//	{
//		int part=4;
//		pw.steps=new ArrayList<Step>();
//		pw.stepsWithOutlier=new ArrayList<Step>();
//                
//        //find first step start
//        int i=pw.getFirstStepIndex();
//		while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//		{
//			if(pw.getY().get(i)>=pw.stepThreshold)
//				break;
//            else
//            	i++;
//		}
//        int app=0;
//        if(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//            app=pw.getY().get(i);
//        else
//            return;
//        for(i=i+1;i<=pw.getLastStepEndIndex() && i<pw.getY().size();)
//        {
//            if(pw.getY().get(i)<app)
//                break;
//            else
//            {
//                app=pw.getY().get(i);
//                i++;
//            }
//        }
//        i--;
//        //looks if it is the real interesting maxima
//        MaximumPos max=new MaximumPos(pw.y.get(i), i);
//        for(int k=i;k<i+pw.distanceAfter && k<=pw.lastStepEndIndex && k<pw.y.size();)
//        {
//            if(pw.y.get(k)>max.getValue())
//                max=new MaximumPos(pw.y.get(k), k);
//            k++;
//        }
//        if(max.getPos()!=i)
//            i=max.getPos();
//        
//        //finds all other steps
//		while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//		{
//			int firstMaxIndex=i;
//			
//			//find stepEquilibrium after the step begin
//			while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//			{
//				if(pw.getY().get(i)<pw.stepEquilibrium)
//					break;
//				else
//					i++;
//			}
//                        
//			//find signal point over stepThreshold
//			while(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//			{
//				if(pw.getY().get(i)>=pw.stepThreshold)
//					break;
//				else
//					i++;
//			}
//			//find maxima
//			if(i<=pw.getLastStepEndIndex() && i<pw.getY().size())
//            {
//                app=pw.getY().get(i);
//                for(i=i+1;i<=pw.getLastStepEndIndex() && i<pw.getY().size();)
//                {
//                    if(pw.getY().get(i)<app)
//                    	break;
//                    else
//                    {
//                        app=pw.getY().get(i);
//                        i++;
//                    }
//                }
//            }
//			
//			//looks if it is the real interesting maxima
//	        max=new MaximumPos(pw.y.get(i), i);
//	        for(int k=i;k<i+pw.distanceAfter && k<=pw.lastStepEndIndex && k<pw.y.size();)
//	        {
//	            if(pw.y.get(k)>max.getValue())
//	                max=new MaximumPos(pw.y.get(k), k);
//	            k++;
//	        }
//	        if(max.getPos()!=i)
//	            i=max.getPos();
//            
//	        //step creation
//			if(i>pw.getLastStepEndIndex())
//				break;
//			else
//			{
//				int overlap=(i-1-firstMaxIndex)/part;
//				if(pw.steps.size()==0)
//					pw.steps.add(new Step(pw,firstMaxIndex, i-1+overlap));
//				else
//					pw.steps.add(new Step(pw,firstMaxIndex-overlap,Integer.min(i-1+overlap, pw.getX().size()-1)));
//			}
//			i--;
//			
//			//avoid loops
//            if(firstMaxIndex==i)
//            	break;
//		}
//	}
	
	/**
	 * It write the step index file in the dataDir folder
	 * @param dataDir directory of the user walk
	 */
	public void printSteps(File dataDir)
	{
		//Write steps in file
		try
		{
			PrintWriter stepWriter=new PrintWriter(new FileOutputStream(dataDir+"\\"+extrapolatedStepsIndexFileName));
			for(Step s : steps)
				stepWriter.write(s.getStartIndexInOriginalFile()+";"+s.getEndIndexInOriginalFile()+"\n");
			stepWriter.close();
		}
		catch (FileNotFoundException e)
		{
			System.out.println("Cannot create file for steps index");
			e.printStackTrace();
		}
		
		try
		{
			PrintWriter stepWriter=new PrintWriter(new FileOutputStream(dataDir+"\\"+extrapolatedStepsWithOutlierIndexFileName));
			for(Step s : stepsWithOutlier)
				stepWriter.write(s.getStartIndexInOriginalFile()+";"+s.getEndIndexInOriginalFile()+"\n");
			stepWriter.close();
		}
		catch (FileNotFoundException e)
		{
			System.out.println("Cannot create file for steps with outlier index");
			e.printStackTrace();
		}
	}
	
	/**
	 * It finds outlier steps, i.e. steps very different from all others found in the walk, and get for free the bestStep, i.e. the step more similar to all others
	 * @param pw walk previously segmented
	 */
	public static void deleteOutlierAndChooseTheBest(PersonWalk pw)
	{
		//prepares data to be DTWed
		int[][] distanceMatrix=new int[pw.getSteps().size()][pw.getSteps().size()];
		int[][][] convertedSteps=new int[pw.getSteps().size()][3][];
		for(int i=0;i<pw.getSteps().size();i++)
		{
			int[] appX=new int[pw.getSteps().get(i).getX().size()];
			int[] appY=new int[pw.getSteps().get(i).getY().size()];
			int[] appZ=new int[pw.getSteps().get(i).getZ().size()];
			for(int j=0;j<appY.length;j++)
			{
				appX[j]=pw.getSteps().get(i).getX().get(j);
				appY[j]=pw.getSteps().get(i).getY().get(j);
				appZ[j]=pw.getSteps().get(i).getZ().get(j);
			}
			convertedSteps[i][0]=appX;
			convertedSteps[i][1]=appY;
			convertedSteps[i][2]=appZ;
		}
		
		//compute DTW distance from steps
		for(int i=0;i<pw.getSteps().size();i++)
		{	
			for(int j=0;j<pw.getSteps().size();j++)
			{
				if(i!=j)
				{
					distanceMatrix[i][j]=(Algorithm.dynamicTimeWarping(convertedSteps[i][0],convertedSteps[j][0])
									+Algorithm.dynamicTimeWarping(convertedSteps[i][1],convertedSteps[j][1])
									+Algorithm.dynamicTimeWarping(convertedSteps[i][2],convertedSteps[j][2]));
				}
				else
					distanceMatrix[i][j]=0;
			}
		}
		
		//compute the avarage distance from step i to all other
		double[] avarageDistanceToAllOther=new double[pw.getSteps().size()];
		for(int i=0;i<avarageDistanceToAllOther.length;i++)
		{
			double app=0;
			for(int j=0;j<pw.getSteps().size();j++)
				app+=distanceMatrix[i][j];
			avarageDistanceToAllOther[i]=app/pw.getSteps().size();
		}
		
		//set the "best" step
		
		//compute medium of medium and variance of steps (+ choosing the best step)
		double distanceThreshold=0;
		int bestStepIndex=0;
		for(int i=0;i<avarageDistanceToAllOther.length;i++)
		{
			distanceThreshold+=avarageDistanceToAllOther[i];
			//check if it is better
			if(avarageDistanceToAllOther[i]<avarageDistanceToAllOther[bestStepIndex])
				bestStepIndex=i;
		}
		distanceThreshold/=avarageDistanceToAllOther.length;
		double variance=0;
		for(double d : avarageDistanceToAllOther)
			variance+=Math.pow((distanceThreshold-d), 2);
		variance/=avarageDistanceToAllOther.length;
		variance=Math.sqrt(variance);
		
		if(bestStepIndex<pw.getSteps().size())
			pw.setBestStep(pw.getSteps().get(bestStepIndex));
		
		//System.out.println("DTW StepThreshold: "+distanceThreshold +" variance: "+variance);
		
		//create "stepsWithOutier" ArrayList before remove steps from "steps"
		pw.setStepsWithOutlier((ArrayList<Step>) pw.getSteps().clone());
		
		//outlier removal
		for(int i=avarageDistanceToAllOther.length-1;i>=0;i--)
		{
			if(avarageDistanceToAllOther[i]>distanceThreshold+variance)
			{
				//remove from "steps"
				pw.getSteps().remove(i);
			}
		}
	}
	
	//It extrapolates all useful data from the probe.
	public static void extrapolateProbe(PersonWalk probe, boolean deleteOutlier)
	{
		/*OLD METHOD
		//separate steps
		separateSteps(probe);
		//delete outlier
		if(deleteOutlier)
			deleteOutlierAndChooseTheBest(probe);
		else
			probe.setStepsWithOutlier((ArrayList<Step>) probe.getSteps().clone());
		*/
		
		probe.stepEquilibrium=(int)(probe.mean-probe.standardDeviation);
        
		//choose threshold
		probe.chooseStepThreshold();

		//separate steps
		separateSteps(probe);

		//delete outlier?
		if(deleteOutlier)
			deleteOutlierAndChooseTheBest(probe);
		else
			probe.setStepsWithOutlier(probe.getSteps());
	}
		
	/**
	 * It loads data previously computed for current user walk
	 * @param filename "data" dir file
	 */
	public void loadSteps(File filename)
	{
		File t=new File(filename+"\\"+stepThresholdFileName);
		if(t.exists())
		{
			try
			{
				stepThreshold=Integer.parseInt(new Scanner(t).nextLine());
			}
			catch (FileNotFoundException e)
			{
				System.out.println("CANNOT FIND \""+stepThresholdFileName+"\" file");
				e.printStackTrace();
			}
		}
		
		t=new File(filename+"\\"+bestStepFileName);
		if(t.exists())
			bestStep=new Step(t);
		else
			System.out.println("CANNOT FIND \""+bestStepFileName+"\" file");
		
		t=new File(filename+"\\"+stepEquilibriumFileName);
		if(t.exists())
		{
			try
			{
				stepEquilibrium=Integer.parseInt(new Scanner(t).nextLine());
			}
			catch (FileNotFoundException e)
			{
				System.out.println("CANNOT FIND \""+stepEquilibriumFileName+"\" file");
				e.printStackTrace();
			}
		}
		
		File f=new File(filename+"\\"+extrapolatedStepsIndexFileName);
		if(f.exists())
		{
			try
			{
				Scanner line=new Scanner(f);
				while(line.hasNextLine())
				{
					String[] app=line.nextLine().split(";");
					steps.add(new Step(this,Integer.parseInt(app[0]),Integer.parseInt(app[1])));
				}
				line.close();
			}
			catch(Exception e)
			{
				System.out.println("CANNOT FIND \""+extrapolatedStepsIndexFileName+"\" file or data are corrupted");
				e.printStackTrace();
			}
		}
		
		f=new File(filename+"\\"+extrapolatedStepsWithOutlierIndexFileName);
		if(f.exists())
		{
			try
			{
				Scanner line=new Scanner(f);
				while(line.hasNextLine())
				{
					String[] app=line.nextLine().split(";");
					stepsWithOutlier.add(new Step(this,Integer.parseInt(app[0]),Integer.parseInt(app[1])));
				}
				line.close();
			}
			catch(Exception e)
			{
				System.out.println("CANNOT FIND \""+extrapolatedStepsWithOutlierIndexFileName+"\" file or data are corrupted");
				e.printStackTrace();
			}
		}
		
		Collections.sort(steps);
		Collections.sort(stepsWithOutlier);
	}
	
	
	//TEST
	
	public static void main(String[] args) throws FileNotFoundException
	{
		//SIMPLE TEST
//		PersonWalk pw=new PersonWalk(new File("PWtry.csv"));
//		System.out.println(pw.getThreshold1step() + " " + pw.firstStepIndex + " " + pw.lastStepEndIndex + " " + pw.mean + " " + pw.standardDeviation);
//		ArrayList<Integer> points=new ArrayList<Integer>();
//		points.add(pw.firstStepIndex);
//		points.add(pw.lastStepEndIndex);
//		SignalShower.genericCurvePlusPoints(pw.getY(), points);
		
//		for(File f : new File("Dataset").listFiles())
//			for(File ff : f.listFiles())
//				for(File fff : ff.listFiles())
//				{
//					PersonWalk pw=new PersonWalk(fff);
//					//System.out.println(pw.getThreshold1step() + " " + pw.firstStepIndex + " " + pw.lastStepEndIndex + " " + pw.mean + " " + pw.standardDeviation);
//					ArrayList<Integer> points=new ArrayList<Integer>();
//					points.add(pw.firstStepIndex);
//					points.add(pw.lastStepEndIndex);
//					//SignalShower.genericCurvePlusPoints(pw.getY(), points);
//				}
	}
	
	//CLASS PARAMETERS
	public static int getSmoothLevel() {return smoothLevel;}
	public static void setSmoothLevel(int smoothLevel) {PersonWalk.smoothLevel = smoothLevel;}
	public static int getSmoothThreshold() {return smoothThreshold;}
	public static void setSmoothThreshold(int smoothThreshold) {PersonWalk.smoothThreshold = smoothThreshold;}
	public int getThresholdLastStepEnd() {return thresholdLastStepEnd;}
	public void setThresholdLastStepEnd(int thresholdLastStepEnd){this.thresholdLastStepEnd = thresholdLastStepEnd;}
	public int getThreshold1step() {return threshold1step;}
	public void setThreshold1step(int threshold1step) {this.threshold1step = threshold1step;}
	//--------------------------------------------------------------------------------
	
	public ArrayList<Step> getSteps() {return steps;}
	public void setSteps(ArrayList<Step> steps) {this.steps = steps;}
	public ArrayList<Step> getStepsWithOutlier() {return stepsWithOutlier;}
	public void setStepsWithOutlier(ArrayList<Step> stepsWithOutlier) {this.stepsWithOutlier = stepsWithOutlier;}
	public int getStepThreshold() {return stepThreshold;}
	public void setStepThreshold(int stepThreshold) {this.stepThreshold = stepThreshold;}
	public Step getBestStep() {return bestStep;}
	public void setBestStep(Step bestStep) {this.bestStep = bestStep;}
	public ArrayList<Integer> getX() {return x;}
	public ArrayList<Integer> getY() {return y;}
	public ArrayList<Integer> getZ() {return z;}
	public ArrayList<Integer> getTime() {return time;}
	public int getFirstStepIndex() {return firstStepIndex;}
	public int getLastStepEndIndex(){return lastStepEndIndex;}
	public File getOriginalFile() {return originalFile;}
	public void setOriginalFile(File originalFile) {this.originalFile = originalFile;}
	public int getStepEquilibrium() {return stepEquilibrium;}
	public void setStepEquilibrium(int stepEquilibrium) {this.stepEquilibrium = stepEquilibrium;}
	public int getStepsPerWalk() {return stepsPerWalk;}
	public void setStepsPerWalk(int stepsPerWalk) {this.stepsPerWalk = stepsPerWalk;}
	public void setX(ArrayList<Integer> x) {this.x = x;}
	public void setY(ArrayList<Integer> y) {this.y = y;}
	public void setZ(ArrayList<Integer> z) {this.z = z;}
	public void setTime(ArrayList<Integer> time) {this.time = time;}		
	
	
	public ArrayList<Double> getArrayDoubleInterpolationX() {
		return arrayDoubleInterpolationX;
	}

	public void setArrayDoubleInterpolationX(ArrayList<Double> arrayDoubleInterpolationX) {
		this.arrayDoubleInterpolationX = arrayDoubleInterpolationX;
	}

	public ArrayList<Double> getArrayDoubleInterpolationY() {
		return arrayDoubleInterpolationY;
	}

	public void setArrayDoubleInterpolationY(ArrayList<Double> arrayDoubleInterpolationY) {
		this.arrayDoubleInterpolationY = arrayDoubleInterpolationY;
	}

	public ArrayList<Double> getArrayDoubleInterpolationZ() {
		return arrayDoubleInterpolationZ;
	}

	public void setArrayDoubleInterpolationZ(ArrayList<Double> arrayDoubleInterpolationZ) {
		this.arrayDoubleInterpolationZ = arrayDoubleInterpolationZ;
	}
	
	
}
