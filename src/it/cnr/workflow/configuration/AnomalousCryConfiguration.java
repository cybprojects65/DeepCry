package it.cnr.workflow.configuration;

public class AnomalousCryConfiguration {


	//################TONE UNIT SEGMENTATION
    public static double maxSilence4ToneUnitSelection = 0.500f;
    public static float minimumAudioLength4ToneUnitSelection = 5; // s
    
    //################ENERGY-PITCH FEATURES SEGMENTATION
    public static float window4EnergyPitchFeatures = 0.1f;
    
    //################ENERGY-PITCH ISLAND SELECTION
    public static double minimumScaledEnergyPitch = 0d;//-0.5d;
	public static int maxTriesToFindValidIslandsAndClusters = 5;
	public static double reductionFactor4Retry = 0.3;
	public static double silenceSecondsToAddBetweenDetectedHighEnergyPitchSegments=0.1;
	
	//################MODULATION SPECTROGRAM
	public static int numberOfMSFeatures = 8;
	public static double maxMSFrequency = 3000;
	
	//(energy> 166 && energy<167 && entropy > 0 && entropy < 1.4) //100% precision
	//################ENTROPY-ENERGY ANOMALY DETECTION
	public static double lowEnergy = 166;
	public static double highEnergy = 167;
	public static double lowEntropy = 0;
	public static double highEntropy = 1.4;
	
	//################DEEP LEARNING
	public static boolean useDLClassification = true;
	public static int nEpochs = 1000;
	public static int reconstructionNumSamples = 16; // fixed
	public static int nhidden = 8;
	public static int nClusters4AnomalyDetection = 5;
	
	//################TEST PARAMETERS
	public static boolean skippreprocessing = false;
	
	//################CRY CLASSIFICATION
	public static double minimumContinuousAnomalousSegment4Classification = 0.0250;//0.0125;//0.200;//0.220 //seconds	
	public static double silenceSecondsToAddBetweenAnomalousCrySamples4Reporting=0.2; 
}
