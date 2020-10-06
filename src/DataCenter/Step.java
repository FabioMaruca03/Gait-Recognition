package DataCenter;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.Scanner;

public class Step implements Comparable<Step>
{
	private ArrayList<Integer> x=new ArrayList<Integer>();
	private ArrayList<Integer> y=new ArrayList<Integer>();
	private ArrayList<Integer> z=new ArrayList<Integer>();
	private ArrayList<Integer> time=new ArrayList<Integer>();
	
	private int startIndexInOriginalFile=-1, endIndexInOriginalFile=-1;

	private int index;
	
	private ArrayList<Double> arrayDoubleInterpolationX;
	private ArrayList<Double> arrayDoubleInterpolationY;
	private ArrayList<Double> arrayDoubleInterpolationZ;
	
	public Step(File f)
	{
		if(f.isFile())
		{
			try
			{
				Scanner in = new Scanner(f);
				while(in.hasNextLine())
				{
					String s[]=in.nextLine().split(";",-1);
					time.add(Integer.parseInt(s[0]));
					x.add(Integer.parseInt(s[1]));
					y.add(Integer.parseInt(s[2]));
					z.add(Integer.parseInt(s[3]));
				}
				in.close();
				if(f.getName().equalsIgnoreCase("BestStep.csv"))
					index=-1;
				else
					index=Integer.parseInt(f.getName().substring(4,f.getName().length()-4));
			}
			catch (FileNotFoundException e)
			{
				System.out.println("Cannot load this step");
				e.printStackTrace();
			}
		}
	}
	
	public Step(ArrayList<Integer> x, ArrayList<Integer> y, ArrayList<Integer> z, int index)
	{
		this.x=x;
		this.y=y;
		this.z=z;
		this.index=index;
	}
	
	public Step(PersonWalk pw, int startStep, int endStep)
	{
		startIndexInOriginalFile=startStep;
		endIndexInOriginalFile=endStep;
		for(int i=startStep;i<=endStep;i++)
		{
			time.add(pw.getTime().get(i));
			x.add(pw.getX().get(i));
			y.add(pw.getY().get(i));
			z.add(pw.getZ().get(i));
		}
	}
	
	public int stepTime()
	{
		return time.get(time.size()-1)-time.get(0);
	}
	
	public String toString()
	{
		String s="";
		int i=0;
		for(;i<y.size()-1;i++)
		{
			s+=time.get(i)+";"+x.get(i)+";"+y.get(i)+";"+z.get(i)+"\n";
		}
		s+=time.get(i)+";"+x.get(i)+";"+y.get(i)+";"+z.get(i);
		return s;
	}

	public ArrayList<Integer> getX() {return x;}
	public void setX(ArrayList<Integer> x) {this.x = x;}
	public ArrayList<Integer> getY() {return y;}
	public void setY(ArrayList<Integer> y) {this.y = y;}
	public ArrayList<Integer> getZ() {return z;}
	public void setZ(ArrayList<Integer> z) {this.z = z;}
	public ArrayList<Integer> getTime() {return time;}
	public void setTime(ArrayList<Integer> time) {this.time = time;}
	public int getStartIndexInOriginalFile() {return startIndexInOriginalFile;}
	public void setStartIndexInOriginalFile(int startIndexInOriginalFile) {this.startIndexInOriginalFile = startIndexInOriginalFile;}
	public int getEndIndexInOriginalFile() {return endIndexInOriginalFile;}
	public void setEndIndexInOriginalFile(int endIndexInOriginalFile) {this.endIndexInOriginalFile = endIndexInOriginalFile;}
	
	@Override
	public int compareTo(Step s)
	{
		return Integer.compare(index, s.index);
	}

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
