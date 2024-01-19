package it.cnr.workflow.coroetal2024.early;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.UUID;

import javax.sound.sampled.AudioFormat;

import it.cnr.clustering.ClassificationManager;
import it.cnr.clustering.DetectionManager;
import it.cnr.clustering.EnergyPitchFilterManager;
import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;
import it.cnr.features.FeatureSet;
import it.cnr.features.IslandDetector;
import it.cnr.models.lstm.DichotomicLSTM;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;
import it.cnr.workflow.utils.UtilsObjects;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class CryDetectorClassifier {

	public MultiKMeans clusterer = new MultiKMeans();
	public DeepLearningManager dlo = new DeepLearningManager();
	public double[][] featureMatrix = null;
	public List<File> extractedWaveFiles = null;
	public List<File> annotationFiles = null;
	public double SNR = 0;
	public File outputFolder;
	public WorkflowConfiguration config;
	
	//################ START OF PARAMETER SECTION ################
	
	//################CRY DETECTION PARAMETERS
	public int maxTriesToFindValidIslandsAndClusters = 5;
	public double reductionFactor4Retry = 0.3;
	public double minimumMSContinuosWindowDetection = 0.5; // seconds
	public int nClassesDetection = 10;
	public int numberOfMSFeatures4Detection = 15;
	
	//################CRY CLASSIFICATION PARAMETERS
	public int nClassesClassification = 40;
	public double entropyEnergyThr = 26;
	public double lowestEntropyEnergyThr = 18;
	public double minimumMSContinuosWindowClassification = 0.3; //seconds
	public int numberOfMSFeatures4Classification = 8;
	public double maxFrequencyForClassification = 3000;
	
	//################ END OF PARAMETER SECTION ################
	
	
	public CryDetectorClassifier() {
		this.config = new WorkflowConfiguration();
		System.out.println("Current configuration:\n" + config.toString());
	}

	public CryDetectorClassifier(WorkflowConfiguration config) {
		this.config = config.clone();
		System.out.println("Current configuration:\n" + config.toString());
	}

	public File detect(File audio) throws Exception {

		// step 1 - separate input file into subfiles based on energy
		System.out.println("#TONE UNIT SEGMENTATION - START#");
		File toneUnitCleanedFile = toneUnitSegmentation(audio);
		System.out.println("#TONE UNIT SEGMENTATION - END#");
		// step 2 - extract features
		System.out.println("#FEATURE EXTRACTION - START#");
		FeatureExtractor featurer = new FeatureExtractor(config);
		// one row for each frame
		double featureMatrix[][] = featurer.extractFeatureMatrix(toneUnitCleanedFile, config.standardiseFeatures);
		System.out.println("#FEATURE EXTRACTION - END#");
		// step 3 - group features in windows

		// step 4 - cluster the features
		DetectionManager clustering = null;
		System.out.println("#CLUSTERING - START#");
		double[] thresholds = UtilsVectorMatrix.initializeVector(featureMatrix[0].length, 1);
		for (int tryn = 1; tryn <= maxTriesToFindValidIslandsAndClusters; tryn++) {
			clustering = new EnergyPitchFilterManager(config, toneUnitCleanedFile);//new DetectionManager(config, toneUnitCleanedFile);
			clustering.detectHighValuedFeatures(featureMatrix, toneUnitCleanedFile.getParentFile(), 1, nClassesDetection, thresholds);

			if (!clustering.highriskclusterfound && tryn == maxTriesToFindValidIslandsAndClusters) {
				System.out.println("A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT!");
				return(null);
			} else if (!clustering.highriskclusterfound) {
				System.out.println(
						"A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT - RETRYING WITH THRESHOLDS REDUCED - RETRY N. "
								+ tryn);
				thresholds = UtilsVectorMatrix.reduceVectorValues(thresholds, reductionFactor4Retry);
				System.out.println("New thresholds: " + Arrays.toString(thresholds));
			} else {
				System.out.println("A HIGH-VALUED CLUSTER WAS FOUND");
				break;
			}

		}
		System.out.println("#CLUSTERING - END#");

		System.out.println("#ISLAND DETECTION - START#");
		IslandDetector islandDetector = new IslandDetector(clustering, config, toneUnitCleanedFile);
		File islandAudio = islandDetector.detectIslands();
		System.out.println("#ISLAND DETECTION - END#");
		
		
		System.out.println("#MODULATION SPECTROGRAM - START#");
		File ms_output = new File(islandAudio.getAbsolutePath().replace(".wav", "_mod_spect.csv"));
		ModulationSpectrogram ms = new ModulationSpectrogram();
		boolean saturateMagnitudeDBs = true;
		boolean addDeltas = true;
		
		double maxMSfrequency = new AudioBits(islandAudio).getAudioFormat().getSampleRate()/2d;
		try {
			ms.calcMS(islandAudio, ms_output, saturateMagnitudeDBs, addDeltas, numberOfMSFeatures4Detection, maxMSfrequency);
		}catch(Exception e) {
			System.out.println("IMPOSSIBLE TO TRACE MODULATION SPECTROGRAM - POSSIBLY BECAUSE CAPTURED AUDIO IS TOO SHORT");
			e.printStackTrace();
			return(null);
		}
		// one row per vector
		// first rows are deltas - next rows double deltas
		//N FEATURES x N VECTORS
		double[][] mod_spec_features = ms.modulationSpectrogram;

		System.out.println("SAVING MS FEATURES");
		File msbinFile = new File(islandAudio.getAbsolutePath().replace(".wav", "ms.bin"));
		FileOutputStream fos = new FileOutputStream(msbinFile);
		ObjectOutputStream oos = new ObjectOutputStream(fos);
		oos.writeObject(mod_spec_features);
		oos.flush();
		oos.close();
		fos.close();
		System.out.println("SAVING MS FEATURES - DONE");
		System.out.println("#MODULATION SPECTROGRAM - END#");

		System.out.println("#RE-CLUSTERING - START#");
		double[][] mod_spec_features_nodelta = UtilsVectorMatrix.subsetRows(mod_spec_features, numberOfMSFeatures4Detection,
				(numberOfMSFeatures4Detection * 2 - 1)); // mod spec without deltas
		WorkflowConfiguration tempConfig = config.clone();
		System.out.println("\tTransposing matrix for clustering");
		//N VECTORS x N FEATURES 
		double[][] mod_spec_features_nodelta_tr = UtilsVectorMatrix.traspose(mod_spec_features_nodelta);
		double[] q3 = UtilsVectorMatrix.columnQ3(mod_spec_features_nodelta_tr);
		System.out.println("\tEnd matrix transposing");

		tempConfig.standardiseFeatures = false;
		DetectionManager clusteringMS = null;
		for (int tryn = 1; tryn <= maxTriesToFindValidIslandsAndClusters; tryn++) {
			clusteringMS = new DetectionManager(tempConfig, islandAudio);
			clusteringMS.detectHighValuedFeatures(mod_spec_features_nodelta_tr, islandAudio.getParentFile(), nClassesDetection, nClassesDetection, q3);
			File labFile = new File(islandAudio.getAbsolutePath().replace(".wav", "modspeclust.lab"));
			AudioBits b = new AudioBits(islandAudio);
			short[] signal = b.getShortVectorAudio();
			int signalLengthMS = signal.length / ModulationSpectrogram.reductionFactor;
			int samplingFreqMS = (int) (b.getAudioFormat().getSampleRate() / ModulationSpectrogram.reductionFactor);
			double windowShiftMS = ModulationSpectrogram.windowShift;
			b.ais.close();
			clusteringMS.toLabCTC(labFile, samplingFreqMS, signalLengthMS, windowShiftMS, minimumMSContinuosWindowDetection);

			if (!clusteringMS.highriskclusterfound && tryn == maxTriesToFindValidIslandsAndClusters) {
				System.out.println("A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT!");
				return(null);
			} else if (!clusteringMS.highriskclusterfound) {
				System.out.println(
						"A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT - RETRYING WITH THRESHOLDS REDUCED - RETRY N. "
								+ tryn);
				q3 = UtilsVectorMatrix.reduceVectorValues(q3, reductionFactor4Retry);
				System.out.println("New thresholds: " + Arrays.toString(q3));
			} else {
				System.out.println("A HIGH-VALUED CLUSTER WAS FOUND");
				break;
			}
		}
		System.out.println("Extracting annotated audio segments");
		File mod_spec_audio = new File(islandAudio.getAbsolutePath().replace(".wav", "_highmodspec.wav"));
		saveNonEmptyAnnotatedSignal(mod_spec_audio, islandAudio, clusteringMS.times, clusteringMS.labels);
		System.out.println("Modulation Spectrogram islands saved to " + mod_spec_audio.getAbsolutePath());
		List<double[]> crySegments = SignalProcessing.getAnnotatedSignalSegments(islandAudio, clusteringMS.times,clusteringMS.labels); 
		File mod_spec_segments = new File(islandAudio.getAbsolutePath().replace(".wav", "_highmodspec.bin"));
		System.out.println("Saving segments to file "+mod_spec_segments.getAbsolutePath());
		saveObject(mod_spec_segments, crySegments);
		System.out.println("#RE-CLUSTERING - END#");		
		
		//TODO: enhance found segments (recall) through LSTM training
		//train on the mod_spec_audio
		//test on the islandAudio
		//record the cry-annotated audio in a new file
		return mod_spec_audio;

	}

	public static void saveNonEmptyAnnotatedSignal(File outputAudiofile, File inputAudioReference, double [] times, String [] labels) throws Exception{
		short[] reducedsignalMS = SignalProcessing.extractAnnotatedSignal(inputAudioReference, times,labels);
		System.out.println("Saving annotated audio segments");
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(reducedsignalMS, outputAudiofile,
			new AudioBits(inputAudioReference).getAudioFormat());
	}
	
	public static void saveObject(File outputFile, Object o) throws Exception {
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outputFile));
		oos.writeObject(o);
		oos.close();
	}

	public static Object loadObject(File inputFile) throws Exception {
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inputFile));
		return (ois.readObject());
	}
	
	public File toneUnitSegmentation(File audio) throws Exception {

		FeatureExtractor extractor = new FeatureExtractor();

		// works with 44100 kHz audio
		outputFolder = extractor.separateFilesBasedOnEnergy(audio, config.maxSilence);
		SNR = extractor.SNR;
		CorpusCleaner.deleteShortAudio(outputFolder, config.minimumAudioLength);
		System.out.println("General SNR " + extractor.SNR);
		File outputCleanedFile = CorpusCleaner.unifyToneUnitsIntoOneFile(outputFolder, audio);

		return outputCleanedFile;

	}




	public void classifyHighModSpecSegments(File audioFiles []) throws Exception {
		
		File mainAudioFolder = audioFiles[0].getParentFile();
		File unifiedAudioFile = new File (mainAudioFolder,"unifiedMS.wav");
		int totalSamples = 0;
		System.out.println("Unifying file to "+unifiedAudioFile.getAbsolutePath());
		List<short[]> signals = new ArrayList<>();
		AudioFormat af = null; 
 		for (File audio:audioFiles) {
 			File audioFolder = new File(audio.getAbsolutePath().replace(".wav", "_processing"));
			File [] allfolderfiles = audioFolder.listFiles();
			File modspecfile = null;
			for (File f:allfolderfiles) {
				if (f.getName().endsWith("_highmodspec.wav")) {
					modspecfile = f;
					break;
				}
			}
			if (modspecfile==null)
				continue;
			
			AudioBits ab = new AudioBits(modspecfile);
			if (af == null) {
				af = new AudioBits(modspecfile).getAudioFormat();
			}
			short[] sig = ab.getShortVectorAudio();
			ab.ais.close();
			totalSamples+=sig.length;
			signals.add(sig);
		}

		short [] totalFile = new short[totalSamples];
		int sam = 0;
		for (short[] sig:signals) {
			for (int i=0;i<sig.length;i++) {
				totalFile[sam] = sig[i];
				sam++;
			}
		}
		
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(totalFile, unifiedAudioFile, af);
 		System.out.println("Unifying file to "+unifiedAudioFile.getAbsolutePath()+" done");
 		File ms_output = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_mod_spect_forClassification.csv"));
 		System.out.println("#MODULATION SPECTROGRAM - START FOR FILE "+unifiedAudioFile.getAbsolutePath());
		ModulationSpectrogram ms = new ModulationSpectrogram();
		boolean saturateMagnitudeDBs = true;
		boolean addDeltas = false;
		ms.calcMS(unifiedAudioFile, ms_output, saturateMagnitudeDBs, addDeltas, numberOfMSFeatures4Classification, maxFrequencyForClassification);
		System.out.println("#MODULATION SPECTROGRAM - END FOR FILE "+unifiedAudioFile.getAbsolutePath());

		System.out.println("#RUNNING CLASSIFICATION");
		double [][] msOverallMatrix = UtilsVectorMatrix.traspose(ms.modulationSpectrogram);
		ClassificationManager classifierMS = new ClassificationManager(config,unifiedAudioFile);
		File labFile = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", ".lab"));
		//classification clustering 
		classifierMS.classifyFeatures(msOverallMatrix, mainAudioFolder, nClassesClassification, nClassesClassification, entropyEnergyThr, lowestEntropyEnergyThr);
		
		int signalLengthMS = totalSamples / ModulationSpectrogram.reductionFactor;
		int samplingFreqMS = (int) (af.getSampleRate() / ModulationSpectrogram.reductionFactor);
		double windowShiftMS = ModulationSpectrogram.windowShift;
		
		classifierMS.toLabCTC(labFile, samplingFreqMS, signalLengthMS, windowShiftMS, minimumMSContinuosWindowClassification);
		
		File mod_spec_audio = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_anomalouscry.wav"));
		saveNonEmptyAnnotatedSignal(mod_spec_audio, unifiedAudioFile, classifierMS.times, classifierMS.labels);
		
		
		System.out.println("#CLASSIFICATION FINISHED");
	}
	
	
  public void classify(File audioFiles []) throws Exception {
		
		File mainAudioFolder = audioFiles[0].getParentFile();
		File unifiedAudioFile = new File (mainAudioFolder,"islands_unified.wav");
		File msFeatureTableObject = new File (unifiedAudioFile.getAbsolutePath().replace(".wav", ".bin"));
		double [][] msOverallMatrix = null;
		AudioFormat af = null; 
		int totalSamples = 0;
		
		if (!msFeatureTableObject.exists()) {
		
		
		System.out.println("Unifying file to "+unifiedAudioFile.getAbsolutePath());
		List<short[]> signals = new ArrayList<>();
		

 		for (File audio:audioFiles) {
 			
 			System.out.println("#ADDING AUDIO FILE "+audio.getName());
 			
 			File audioFolder = new File(audio.getAbsolutePath().replace(".wav", "_processing"));
			File [] allfolderfiles = audioFolder.listFiles();
			File modspecfile = null;
			for (File f:allfolderfiles) {
				if (f.getName().endsWith("_islands.wav")) {
					modspecfile = f;
					break;
				}
			}
			if (modspecfile==null)
				continue;
			
			AudioBits ab = new AudioBits(modspecfile);
			if (af == null) {
				af = new AudioBits(modspecfile).getAudioFormat();
			}
			short[] sig = ab.getShortVectorAudio();
			ab.ais.close();
			totalSamples+=sig.length;
			signals.add(sig);
		}

		short [] totalFile = new short[totalSamples];
		int sam = 0;
		for (short[] sig:signals) {
			for (int i=0;i<sig.length;i++) {
				totalFile[sam] = sig[i];
				sam++;
			}
		}
		
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(totalFile, unifiedAudioFile, af);
 		System.out.println("Unifying file to "+unifiedAudioFile.getAbsolutePath()+" done");
 		File ms_output = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_mod_spect_forClassification.csv"));
 		System.out.println("#MODULATION SPECTROGRAM - START FOR FILE "+unifiedAudioFile.getAbsolutePath());
		ModulationSpectrogram ms = new ModulationSpectrogram();
		boolean saturateMagnitudeDBs = true;
		boolean addDeltas = false;
		ms.calcMS(unifiedAudioFile, ms_output, saturateMagnitudeDBs, addDeltas, numberOfMSFeatures4Classification, maxFrequencyForClassification);
		System.out.println("#MODULATION SPECTROGRAM - END FOR FILE "+unifiedAudioFile.getAbsolutePath());
		msOverallMatrix = UtilsVectorMatrix.traspose(ms.modulationSpectrogram);
			saveObject(msFeatureTableObject, msOverallMatrix);
		}else {
			msOverallMatrix = (double [][]) loadObject(msFeatureTableObject);
			AudioBits ab = new AudioBits(unifiedAudioFile);
			af = ab.getAudioFormat();
			short[] sig = ab.getShortVectorAudio();
			ab.ais.close();
			totalSamples=sig.length;
		}
		
		System.out.println("#RUNNING CLASSIFICATION");
		
		ClassificationManager classifierMS = new ClassificationManager(config,unifiedAudioFile);
		File labFile = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", ".lab"));
		//classification clustering 
		classifierMS.classifyFeatures(msOverallMatrix, mainAudioFolder, nClassesClassification, nClassesClassification, entropyEnergyThr, lowestEntropyEnergyThr);
		
		int signalLengthMS = totalSamples / ModulationSpectrogram.reductionFactor;
		int samplingFreqMS = (int) (af.getSampleRate() / ModulationSpectrogram.reductionFactor);
		double windowShiftMS = ModulationSpectrogram.windowShift;
		
		classifierMS.toLabCTC(labFile, samplingFreqMS, signalLengthMS, windowShiftMS, minimumMSContinuosWindowClassification);
		
		File mod_spec_audio = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_anomalouscry.wav"));
		saveNonEmptyAnnotatedSignal(mod_spec_audio, unifiedAudioFile, classifierMS.times, classifierMS.labels);
		
		FeatureSet fs = FeatureSet.fromModulationSpectrogram(msOverallMatrix, classifierMS.times, classifierMS.labels, windowShiftMS);
		File featureSetFile = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_featureset.bin"));
		saveObject(featureSetFile, fs);
		
		System.out.println("#CLASSIFICATION FINISHED");
		
		System.out.println("#EXTENDING CLASSIFICATION");
		System.out.println("Extending the times of the found MS regions");
		IslandDetector id = new IslandDetector(config, unifiedAudioFile);
		List<double []> times = id.extendIntervalsThroughEnergyIslands(fs.timeIntervals);
		double tlast = 0;
		List<Double> timesN = new ArrayList<>();
		List<String> labelsN = new ArrayList<>();
		timesN.add(0d);
		labelsN.add("");
		  
		for (double [] tint : times) {
			  
			  double t0 = tint[0];
			  double t1 = tint[1];
			  String l = "A";
			  if (t0!=tlast) {
				  timesN.add(t0);
				  labelsN.add(l);
				  timesN.add(t1);
				  labelsN.add("");
			  }
		  }
		  
		  double[] timesDA = new double[timesN.size()];
		  for (Integer h = 0; h <timesDA.length;h++) {
			  timesDA[h] = timesN.get(h);
		  }

		  String[] labelsDA = new String[labelsN.size()];

		  for (Integer h = 0; h <timesDA.length;h++) {
			  labelsDA[h] = labelsN.get(h);
		  }
		  
		  File labeloutput = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_islcl.lab"));
		  CorpusCleaner.vector2LabFile(labelsDA, timesDA, labeloutput);
		
		  System.out.println("#EXTENDING CLASSIFICATION FINISHED");
		  
		  System.out.println("#DONE");
	}

  public void extendClassification(File msFeatureTableObjectFile, File featureSetObjectFile, File unifiedAudioFile) throws Exception {
	  
	  System.out.println("LOADING MODULATION SPECTROGRAM MATRIX");
	  double [][] msOverallMatrix = (double [][]) loadObject(msFeatureTableObjectFile);
	  
	  System.out.println("LOADING CHARACTERISED CRY-RELATED FEATURE SET MATRIX");
	  FeatureSet fs = (FeatureSet)  loadObject(featureSetObjectFile);
	  
	  System.out.println("#COMPLETE MODULATION SPECTROGRAM - START FOR FILE "+unifiedAudioFile.getAbsolutePath());
	  File ms_output_complete = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_complete_mod_spect_forClassification.csv"));
	  File ms_output_completebin = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_complete_mod_spect_forClassification.bin"));
	  double[][] mod_spec_complete = null;
	  if (!ms_output_completebin.exists()) {
		  ModulationSpectrogram ms_complete = new ModulationSpectrogram();
		  double maxMSfrequency = new AudioBits(unifiedAudioFile).getAudioFormat().getSampleRate()/2d;
		  ms_complete.calcMS(unifiedAudioFile, ms_output_complete, true, false, numberOfMSFeatures4Detection, maxMSfrequency);
		  mod_spec_complete = UtilsVectorMatrix.traspose(ms_complete.modulationSpectrogram);
		  saveObject(ms_output_completebin, mod_spec_complete);
	  }else {
		  mod_spec_complete = (double[][]) loadObject(ms_output_completebin);
	  }
	  System.out.println("#COMPLETE MODULATION SPECTROGRAM - END FOR FILE "+unifiedAudioFile.getAbsolutePath());
		
		double[] msq3 = UtilsVectorMatrix.columnQ3(mod_spec_complete);
	  File tempFolder = new File(unifiedAudioFile.getParentFile() ,"temp"+UUID.randomUUID());

	  System.out.println("creating temporary folder "+tempFolder.getAbsolutePath());
	  tempFolder.mkdir();
	  File trainingFolder = new File(tempFolder,"train");
	  File testFolder = new File(tempFolder,"test");
	  trainingFolder.mkdir();
	  testFolder.mkdir();
	  File trainingFeaturesFolder = new File(trainingFolder,"features");	  
	  trainingFeaturesFolder.mkdir();
	  File trainingLabelsFolder = new File(trainingFolder,"labels");
	  trainingLabelsFolder.mkdir();
	  File testFeaturesFolder = new File(testFolder,"features");
	  testFeaturesFolder.mkdir();
	  File testLabelsFolder = new File(testFolder,"labels");
	  testLabelsFolder.mkdir();
	  int counter = 0;
	  counter = UtilsObjects.saveTS4DL(fs.featuresPerTimeInterval,1,counter,trainingFeaturesFolder,testFeaturesFolder, trainingLabelsFolder,  testLabelsFolder);
	  
	  DichotomicLSTM lstm = new DichotomicLSTM();
		lstm.init(tempFolder, config.nhidden, config.nClasses, config.minibatch, config.nEpochs);
		//prepare the data 
		lstm.train();
		System.out.println("LSTM accuracy: "+lstm.accuracy);
		System.out.println("LSTM evaluation: "+lstm.evaluation);
		
	  DichotomicLSTM lstmclassification = new DichotomicLSTM();
	  lstmclassification.init(tempFolder, config.nhidden, config.nClasses, config.minibatch, config.nEpochs);
	  
	  System.out.println("Classification of the training set");
	  int i = 0;
	  for (double ff[][]:fs.featuresPerTimeInterval) {
		  int c = lstmclassification.classify(ff);
		  System.out.println("Set n. "+i+" classified as "+c);
	  }
	  
	  int itShift = 1;
	  int nelementsPerTimeSeries = ModulationSpectrogram.getMSIndex(minimumMSContinuosWindowClassification, ModulationSpectrogram.windowShift)+1;
	  int it0 = 0;
	  int it1 = it0+nelementsPerTimeSeries;

	  double classifThreshold = 0.70;
	  
	  int max_it = ModulationSpectrogram.getMSIndex(20, ModulationSpectrogram.windowShift)+1;
	  max_it = msOverallMatrix.length;
	  int ncol = msOverallMatrix[0].length;
	  int ncolComplete = mod_spec_complete[0].length;
	  
	  //LinkedHashSet<Double> times = new LinkedHashSet<Double>();
	  //LinkedHashSet<String> labels = new LinkedHashSet<String>();
	  //List<Double> times = new ArrayList<Double>();
	  //List<String> labels = new ArrayList<String>();
	  HashMap<Integer, Double> times = new LinkedHashMap<>();
	  HashMap<Integer, String> labels = new LinkedHashMap<>();
	  
	  boolean isRecording = false;
	  int sid = 0;
	  while (it1<max_it) {
		  double t0 = it0*ModulationSpectrogram.windowShift;
		  double t1 = it1*ModulationSpectrogram.windowShift;
		  
		  double [][] submatrix = new double [it1-it0][ncol];
		  for (int r=it0;r<it1;r++) {
			  for (int c=0;c<ncol;c++) {
				  submatrix[r-it0][c] = msOverallMatrix[r] [c];
			  }
		  }
		  
		  int cl = lstmclassification.classify(submatrix);
		  if (lstmclassification.class1Score>classifThreshold ) {
			  
			  double [][] submatrixComplete = new double [it1-it0][ncolComplete];
			  for (int r=it0;r<it1;r++) {
				  for (int c=0;c<ncolComplete;c++) {
					  submatrixComplete[r-it0][c] = mod_spec_complete[r] [c];
				  }
			  }
			  
			  double cm [] = UtilsVectorMatrix.columnMeans(submatrixComplete);
			  int cmhigh = 0;
			  for (int v=0;v<cm.length;v++) {
				  if (cm[v]>=msq3[v])
						  cmhigh++;
			  }
			  boolean highMS = false;
			  double ratio = UtilsMath.roundDecimal(((double) cmhigh/ (double) cm.length) ,2);
					  
			  if (cmhigh>0.5*cm.length)
				  highMS = true;
			  
		  //System.out.println("Set n. "+i+" classified as "+c);
			  //System.out.println("["+Utils.roundDecimal(t0,2)+";"+Utils.roundDecimal(t1,2)+"] "+1 +" ("+lstmclassification.class1Score+"; HMS:"+highMS+ " $ "+ratio+")");
			  		  
			  if (!isRecording && highMS) {
				  System.out.println("["+UtilsMath.roundDecimal(t0,2)+";"+UtilsMath.roundDecimal(t1,2)+"] "+1 +" ("+lstmclassification.class1Score+"; HMS:"+highMS+ " $ "+ratio+")");
				  isRecording = true;
				  times.put(sid,t0);
				  labels.put(sid, "A");
				  //times.add(t0);
				  //labels.add("Anomalous DP");
				  sid++;
			  }
			  
			  it0 = it1;
			  //it0 = it0+itShift;
			  it1 = it0+nelementsPerTimeSeries;
			  
		  }else {
			  if (isRecording) {
				  System.out.println("t_final = "+t1);
				  isRecording = false;
				  //times.add(t1);
				  //labels.add("");
				  times.put(sid,t1);
				  labels.put(sid, "");
				  sid++;
			  }
			  it0 = it0+itShift;
			  it1 = it0+nelementsPerTimeSeries;
		  }
		  
		  
	  }
	
	  System.out.println("Producing times and labels");
	  
	  double[] timesDA = new double[times.size()];
	  //int h = 0;
	  //for (Double d:times) {
	  for (Integer h = 0; h <timesDA.length;h++) {
		  timesDA[h] = times.get(h);
	  }

	  String[] labelsDA = new String[labels.size()];

	  for (Integer h = 0; h <timesDA.length;h++) {
		  labelsDA[h] = labels.get(h);
	  }
	  
	  File labeloutput = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_dl.lab"));
	  CorpusCleaner.vector2LabFile(labelsDA, timesDA, labeloutput);
	  
	  System.out.println("DONE.");
	  
  }

 public void extendClassificationIslands(File msFeatureTableObjectFile, File featureSetObjectFile, File unifiedAudioFile) throws Exception {
	  
	  System.out.println("LOADING MODULATION SPECTROGRAM MATRIX");
	  double [][] msOverallMatrix = (double [][]) loadObject(msFeatureTableObjectFile);
	  
	  System.out.println("LOADING CHARACTERISED CRY-RELATED FEATURE SET MATRIX");
	  FeatureSet fs = (FeatureSet)  loadObject(featureSetObjectFile);
	  
	  System.out.println("Extending the times of the found MS regions");
	  IslandDetector id = new IslandDetector(config, unifiedAudioFile);
	  List<double []> times = id.extendIntervalsThroughEnergyIslands(fs.timeIntervals);
	  
	  double tlast = 0;
	  
	  List<Double> timesN = new ArrayList<>();
	  List<String> labelsN = new ArrayList<>();
	  timesN.add(0d);
	  labelsN.add("");
	  
	  for (double [] tint : times) {
		  
		  double t0 = tint[0];
		  double t1 = tint[1];
		  String l = "A";
		  if (t0!=tlast) {
			  timesN.add(t0);
			  labelsN.add(l);
			  timesN.add(t1);
			  labelsN.add("");
		  }
	  }
	  
	  double[] timesDA = new double[timesN.size()];
	  for (Integer h = 0; h <timesDA.length;h++) {
		  timesDA[h] = timesN.get(h);
	  }

	  String[] labelsDA = new String[labelsN.size()];

	  for (Integer h = 0; h <timesDA.length;h++) {
		  labelsDA[h] = labelsN.get(h);
	  }
	  
	  File labeloutput = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_islcl.lab"));
	  CorpusCleaner.vector2LabFile(labelsDA, timesDA, labeloutput);
	  
	  System.out.println("DONE.");
	  
  }

}
