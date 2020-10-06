package DataCenter;

import java.io.File;
import java.util.ArrayList;

import ImageAnalysis.Constants;
import ImageAnalysis.ImageCreator;
import ImageAnalysis.UserWalk;
import Utilities.Algorithm;

//import Recognizer.GraphicRecognizer;
//import Utilities.DatasetConverter;

public class PeopleDataset
{
	ArrayList<PersonData> users=new ArrayList<PersonData>();
	boolean errorDuringLoad=false;

	public PeopleDataset(File folder)
	{
//		if (folder != null) {
//			System.out.print(folder.getName());
//			if (folder.isDirectory()) {
//				System.out.print("YES");
//			} else {
//				System.out.print("NO");
//			}
//		}
		
		
		if(folder.isDirectory())
		{
			try
			{
				File people[]=folder.listFiles();
				for(File p:people)
				{
					if(p != null) {
						System.out.println(p.getName());
					}
					add(p);
				}
			}
			catch(Exception e)
			{
				e.printStackTrace();
				//GraphicRecognizer.showErrorInStepSegmentation();
				return;
			}
		}
		else
		{
			System.out.println("Is not a Dataset folder; System cannot load data");
			errorDuringLoad=true;
		}
		
		
	}

	public PeopleDataset()
	{

	}

	public PeopleDataset(ArrayList<PersonData> dataset)
	{
		users=dataset;
	}

	public void add(PersonData pd)
	{
		users.add(pd);
	}

	public void add(File filename)
	{
		PersonData p=new PersonData(filename);
		if(!p.errorDuringLoad)
			users.add(p);
		else
			errorDuringLoad=true;
	}

	public String toString()
	{
		String dataset="";
		for(PersonData pd:users)
			dataset+=pd;
		return dataset;
	}

	public ArrayList<PersonData> getUsers() {return users;}
	public boolean getErrorDuringLoad() {return errorDuringLoad;}

	//TEST
	public static void main(String[] args)
	{
		
		//Lettura dati dal file CVS e immissione nel'oggetto peopleDataset ()
		
		PeopleDataset peopleDataset = new PeopleDataset(new File("DatasetChinaComplete-Segmented"));

		for (PersonData personData:peopleDataset.getUsers()) {
			System.out.println(personData.getPersonWalks().size());
			for (PersonWalk personWalk:personData.getPersonWalks()) {
				System.out.println("personWalk.ToString: " + personWalk.toString());
				personWalk.setArrayDoubleInterpolationX(Algorithm.linearInterpolation(personWalk.getX(), personWalk.getTime(), 1, null));
				personWalk.setArrayDoubleInterpolationY(Algorithm.linearInterpolation(personWalk.getY(), personWalk.getTime(), 1, null));
				personWalk.setArrayDoubleInterpolationZ(Algorithm.linearInterpolation(personWalk.getZ(), personWalk.getTime(), 1, null));
				
				System.out.println("getArrayDoubleInterpolationX " + personWalk.getArrayDoubleInterpolationX());
				System.out.println("getArrayDoubleInterpolationY " + personWalk.getArrayDoubleInterpolationY());
				System.out.println("getArrayDoubleInterpolationZ " + personWalk.getArrayDoubleInterpolationZ());
				System.out.println("steps " + personWalk.getSteps().size());
				for (Step step:personWalk.getSteps()) {
					//dt = 1, future = 0
					
					step.setArrayDoubleInterpolationX(Algorithm.linearInterpolation(step.getTime(), step.getX(), 1, null));
					step.setArrayDoubleInterpolationY(Algorithm.linearInterpolation(step.getTime(), step.getY(), 1, null));
					step.setArrayDoubleInterpolationZ(Algorithm.linearInterpolation(step.getTime(), step.getZ(), 1, null));
					ArrayList<Double> arrayDoubleX = Algorithm.interpolatedNoiseRemoval(step.getTime(), step.getX(), 1, 0);
					ArrayList<Double> arrayDoubleY = Algorithm.interpolatedNoiseRemoval(step.getTime(), step.getY(), 1, 0);
					ArrayList<Double> arrayDoubleZ = Algorithm.interpolatedNoiseRemoval(step.getTime(), step.getZ(), 1, 0);
				}
			}
		}
		
		System.out.println("people in dataset: " + peopleDataset.getUsers().size());

		
		for (PersonData personData:peopleDataset.getUsers()) {
			int allStepsLength = 0;
			int totalSteps = 0;
			for (PersonWalk walk : personData.getPersonWalks()) {
				totalSteps += walk.getSteps().size();
				
				for (Step step : walk.getSteps()) {
					System.out.println("Step: " + step);
					allStepsLength += (step.getEndIndexInOriginalFile() + 1);
				}
			}
			System.out.println("allStepsLength " + allStepsLength);
			
			// calcola la lunghezza media di ogni passo
			double meanStepLength = (double)allStepsLength / totalSteps;
			
			
			personData.allStepsLength = allStepsLength;
			personData.totalSteps = totalSteps;
			personData.meanStepLength = meanStepLength;
		}
		
		for (PersonData personData:peopleDataset.getUsers()) {
			System.out.println("personData: " + personData.getName());
			System.out.println("personData.allStepsLength " + personData.allStepsLength);
			System.out.println("personData.totalSteps" + personData.totalSteps);
			System.out.println("personData.meanStepLength " + personData.meanStepLength + "\n");
		}



//		ImageCreator imageCreator = new ImageCreator();
		
//		for (PersonData personData:peopleDataset.getUsers()) {
//			
////			int walkId = 0;
////			for (PersonWalk walk : personData.getPersonWalks()) {
////				System.out.println(personData.toString());
////				// crea l'immagine
////				System.out.println(" WIDTH: " + personData.meanStepLength);
////				imageCreator.createImage2(Math.ceil(personData.meanStepLength), Constants.datasetSource + personData.getName() +File.pathSeparator+"walk" + walkId, walk, walkId);
////				walkId++;
////			}
//		}
		System.out.println("FINE");
	}
}
