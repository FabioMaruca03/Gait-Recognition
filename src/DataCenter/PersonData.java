package DataCenter;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

import static AppClasses.Constants.*;

public class PersonData
{
	
	int allStepsLength = 0;
	int totalSteps = 0;
	double meanStepLength=0.0;
	
	private ArrayList<PersonWalk> personWalks = new ArrayList<PersonWalk>();
	private String name;
	boolean errorDuringLoad=false;
	
	public PersonData(ArrayList<PersonWalk> pw, String name)
	{
		this.name=name;
		personWalks=pw;
	}
	
	public PersonData(File filename)
	{
		name=filename.getName();
		if(filename.isDirectory())
		{
			File files[]=filename.listFiles();
			for(File f:files)
			{
				
				if(f.isDirectory())
				{
					PersonWalk walk=null;
					for(File app : f.listFiles())
					{
						if(app.isFile())
						{
							walk=new PersonWalk(app);
						}
					}
					if(walk!=null)
					{
						boolean foundData=false;
						for(File app : f.listFiles())
						{
							if(app.isDirectory())
							{
								File st=new File(app.getAbsolutePath()+"\\"+stepThresholdFileName);
								File se=new File(app.getAbsolutePath()+"\\"+stepEquilibriumFileName);
								File bs=new File(app.getAbsolutePath()+"\\"+bestStepFileName);
								File sd=new File(app.getAbsolutePath()+"\\"+extrapolatedStepsIndexFileName);
								File swod=new File(app.getAbsolutePath()+"\\"+extrapolatedStepsWithOutlierIndexFileName);
								if(st.exists() && bs.exists() && sd.exists() && swod.exists() && se.exists())
								{
									try
									{
										walk.setStepThreshold(Integer.parseInt(new Scanner(st).nextLine()));
									}
									catch (Exception e) {e.printStackTrace();}
									try
									{
										walk.setStepEquilibrium(Integer.parseInt(new Scanner(se).nextLine()));
									}
									catch (Exception e) {e.printStackTrace();}
									
									walk.setBestStep(new Step(bs));
									walk.loadSteps(app);
								}
								else
									walk.createSteps(f);
								foundData=true;
							}
						}
						if(!foundData)
							walk.createSteps(f);
						personWalks.add(walk);
					}
				}
			}
		}
		else
		{
			System.out.println(name + " Is not a folder; System cannot load data");
			errorDuringLoad=true;
		}
		
	}
	
	private void loadWalks(File f)
	{
		File files[]=f.listFiles();
		for(File w:files)
			add(w);
	}

	public void add(PersonWalk walk)
	{
		personWalks.add(walk);
	}
	
	public void add(File filename)
	{
		personWalks.add(new PersonWalk(filename));
	}
	
	public String toString()
	{
		String person=name+"\n";
		for(int i=0;i<personWalks.size();i++)
			person+="Walk #"+i+1+":\n"+personWalks.get(i);
		
		return person;
	}

	public ArrayList<PersonWalk> getPersonWalks() {return personWalks;}
	public String getName() {return name;}
}
