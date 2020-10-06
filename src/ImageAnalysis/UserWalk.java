package ImageAnalysis;

import java.util.List;

import org.opencv.core.Mat;

public class UserWalk {

	private int user;
	private int walk;
	private List<Integer> stepsStartIndex;
	private List<Integer> stepsEndIndex;
	private List<Step> stepsList;
	private List<Mat> histograms;

	public int getUser() {return user;}
	public void setUser(int user) {this.user = user;}

	public int getWalk() {return walk;}
	public void setWalk(int walk) {this.walk = walk;}
	
	public List<Integer> getStepsStartIndex() {return stepsStartIndex;}
	public void setStepsStartIndex(List<Integer> startIndex) {this.stepsStartIndex = startIndex;}
	
	public List<Integer> getEndIndex() {return stepsEndIndex;}
	public void setEndIndex(List<Integer> endIndex) {this.stepsEndIndex = endIndex;}
	
	public List<Step> getStepsList() {return stepsList;}
	public void setStepsList(List<Step> stepList) {this.stepsList = stepList;}
	
	public List<Mat> getHistograms() { return histograms; }
	public void setHistograms(List<Mat> histograms) { this.histograms = histograms; }

}
