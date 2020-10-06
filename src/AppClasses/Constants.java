package AppClasses;

import java.io.File;

/**
 * Class containing constants
 * NOTE: When you create a new recognition method you have to add the corresponding "%ModeTestFile", "%ModeTestTableFile", "%TestMode", "%AlgoName" constants.
 * Moreover, you have to add the "%TestMode" variable to the modeArray and the "%AlgoName" to the algoNames array.
 * % represents the name of your new recognition method
 * @author aless
 *
 */
public class Constants
{
	public static final String extrapolatedDataFolderName="data";
	public static final String extrapolatedStepsIndexFileName="stepsIndex.txt";
	public static final String extrapolatedStepsWithOutlierIndexFileName="stepsWithOutlierIndex.txt";
	
	public static final String KDatasetFolder="KDataset";
	
	public static final String walks="walks";
	public static final String stepThresholdFileName="StepsThreshold.txt";
	public static final String bestStepFileName="BestStep.csv";
	public static final String stepEquilibriumFileName="stepEquilibrium.txt";
	
	public static final String EERFolder="EER";
	public static final String WEER=EERFolder+"\\Walk_EER.txt";
	public static final String BSEER=EERFolder+"\\BestStep_EER.txt";
	public static final String BSVSAEER=EERFolder+"\\BestStepVSAll_EER.txt";
	public static final String ASVSAEER=EERFolder+"\\AllStepsVSAll_EER.txt";
	public static final String SSWEER=EERFolder+"\\StepsSW_EER.txt";
	
	// SINGLE TEST FOLDER / SINGLE TEST FILE NAME
	public static final String testDatasetResultFolder="TestDatasetResult";
	public static final File walkModeTestFile=new File(testDatasetResultFolder+"\\WalkModeTestDatasetResult.txt");
	public static final File bestStepModeTestFile=new File(testDatasetResultFolder+"\\BestStepModeTestDatasetResult.txt");
	public static final File bestStepVSAllModeTestFile=new File(testDatasetResultFolder+"\\BestStepVSAllModeTestDatasetResult.txt");
	public static final File allStepsVSAllModeTestFile=new File(testDatasetResultFolder+"\\AllStepsVSAllModeTestDatasetResult.txt");
	public static final File stepSlidingWindowModeTestFile=new File(testDatasetResultFolder+"\\StepSlidingWindowModeTestDatasetResult.txt");
	public static final File magnitudeWalkModeTestFile=new File(testDatasetResultFolder+"\\MagnitudeWalkModeTestDatasetResult.txt");
	public static final File HistogramDistanceWalkModeTestFile=new File(testDatasetResultFolder+"\\HistogramDistanceWalkModeTestDatasetResult.txt");
	public static final File WalkSkylineModeTestFile=new File(testDatasetResultFolder+"\\WalkSkylineModeTestDatasetResult.txt");
	
	//TEST FOLDER / TEST TABLE FILE NAMES
	public static final String testTableResultFolder="TestTableResult";
	public static final File walkModeTestTableFile=new File(testTableResultFolder+"\\WalkModeTestDatasetResult.csv");
	public static final File bestStepModeTestTableFile=new File(testTableResultFolder+"\\BestStepModeTestDatasetResult.csv");
	public static final File bestStepVSAllModeTestTableFile=new File(testTableResultFolder+"\\BestStepVSAllModeTestDatasetResult.csv");
	public static final File allStepsVSAllModeTestTableFile=new File(testTableResultFolder+"\\AllStepsVSAllModeTestDatasetResult.csv");
	public static final File stepSlidingWindowModeTestTableFile=new File(testTableResultFolder+"\\StepSlidingWindowModeTestDatasetResult.csv");
	public static final File magnitudeWalkModeTestTableFile=new File(testTableResultFolder+"\\MagnitudeWalkModeTestDatasetResult.csv");
	public static final File HistogramDistanceWalkModeTestTableFile=new File(testTableResultFolder+"\\HistogramDistanceWalkModeTestDatasetResult.csv");
	public static final File WalkSkylineModeTestTableFile=new File(testTableResultFolder+"\\WalkSkylineModeTestDatasetResult.csv");
	
	
	public static final String defaultProbeLocation="\\Probe";
	public static final String defaultDatasetLocation="\\Dataset";
	
	
	public static final String bestWeightFolder="BestWeight";
	public static final String bestWalkModeWeightsFile=bestWeightFolder+"\\Best-WalkMode-Weights.txt";
	public static final String bestBestStepModeWeightsFile=bestWeightFolder+"\\Best-BestStepMode-Weights.txt";
	public static final String bestBestStepVSAllModeWeightsFile=bestWeightFolder+"\\Best-BestStepVSAllMode-Weights.txt";
	public static final String bestAllStepsVSAllModeWeightsFile=bestWeightFolder+"\\Best-AllStepsVSAllMode-Weights.txt";
	public static final String bestStepSlidingWindowModeWeightsFile=bestWeightFolder+"\\Best-StepSlidingWindowMode-Weights.txt";
	
	public static final String bestPhiFolder="BestPhi";
	public static final String bestPhiWalkModeFile=bestPhiFolder+"\\bestPhiWalkMode.txt";
	public static final String bestPhiBestStepModeFile=bestPhiFolder+"\\bestPhiBestStepMode.txt";
	public static final String bestPhiBestStepVSAllModeFile=bestPhiFolder+"\\bestPhiBestStepVSAllMode.txt";
	public static final String bestPhiAllStepsVSAllModeFile=bestPhiFolder+"\\bestPhiAllStepsVSAllMode.txt";
	public static final String bestPhiStepSlidingWindowModeFile=bestPhiFolder+"\\bestPhiStepSlidingWindowMode.txt";
	
	public static final String openSetThresholdFolder="OpenSetThreshold";
	public static final String walkModeOSTFile=openSetThresholdFolder+"\\walkModeOST.txt";
	public static final String bestStepModeOSTFile=openSetThresholdFolder+"\\bestStepModeOST.txt";
	public static final String bestStepVSAllModeOSTFile=openSetThresholdFolder+"\\bestStepVSAllModeOST.txt";
	public static final String allStepsVSAllModeOSTFile=openSetThresholdFolder+"\\allStepsVSAllModeOST.txt";
	public static final String stepsSlidingWindowModeOSTFile=openSetThresholdFolder+"\\stepsSlidingWindowModeOST.txt";
	
	//TEST MODE - ALGORITHM
	public static final int walkTestMode=0;
	public static final int bestStepTestMode=1;
	public static final int bestStepVSAllTestMode=2;
	public static final int allStepsVSAllTestMode=3;
	public static final int stepSlidingWindowTestMode=4;
	public static final int magnitudeWalkTestMode=5;
	public static final int histogramDistanceWalkTestMode=6;
	public static final int walkSkyLineTestMode=7;
	public static final int testAllMode=999; //used for setup methods
	
	//TEST MODE - ALGORITHM NAMES
	public static final String walkAlgoName="Walk";
	public static final String bestStepAlgoName="Best Step";
	public static final String bestStepVSAllAlgoName="Best Step vs. All";
	public static final String allStepsVSAllAlgoName="All Steps vs. All";
	public static final String stepSlidingWindowAlgoName="Step Sliding Window";
	public static final String magnitudeWalkAlgoName="Magnitude Walk";
	public static final String histogramDistanceWalkAlgoName="Histogram Similarity";
	public static final String walkSkylineAlgoName="Walk Skyline Similarity";
	
	// ARRAY with all mode/algorithm names
	public static final int[] modeArray={walkTestMode,
										 bestStepTestMode,
										 bestStepVSAllTestMode,
										 allStepsVSAllTestMode,
										 stepSlidingWindowTestMode,
										 magnitudeWalkTestMode,
										 histogramDistanceWalkTestMode,
										 walkSkyLineTestMode};
	
	public static final String[] algoNames={walkAlgoName,
											bestStepAlgoName,
											bestStepVSAllAlgoName,
											allStepsVSAllAlgoName,
											stepSlidingWindowAlgoName,
											magnitudeWalkAlgoName,
											histogramDistanceWalkAlgoName,
											walkSkylineAlgoName};
	
	public static final int identificationCSModality=0;
	public static final int identificationOSModality=1;
	public static final int verificationModality=2;
	
	public static final String IDCSMCode="IDCS";
	public static final String IDOSMCode="IDOS";
	public static final String VERMCode="VER";
	
	public static final int minimumStepSize=40;
	
	//It is used in order to decide the approximation of the threshold for the test in verification (the threshold goes from [0,1] with an increment of 1/defaultVerifyGranularity)
	public static final int defaultVerifyGranularity=1000;
	
	
	public static final int identificationCSThreshold=Integer.MAX_VALUE;
}
