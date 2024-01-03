package it.cnr.workflow.coroetal2024;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import javax.sound.sampled.AudioFormat;

import org.apache.commons.io.FileUtils;

import it.cnr.clustering.ClassificationManager;
import it.cnr.clustering.ClusteringManager;
import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;
import it.cnr.features.IslandDetector;
import it.cnr.features.Utils;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.speech.filters.GreenwoodFilterBank;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.speech.filters.SignalConverter;
import it.cnr.workflow.configuration.Configuration;

public class CryDetectorClassifier {

	public MultiKMeans clusterer = new MultiKMeans();
	public DeepLearningManager dlo = new DeepLearningManager();
	public double[][] featureMatrix = null;
	public List<File> extractedWaveFiles = null;
	public List<File> annotationFiles = null;
	public double SNR = 0;
	public File outputFolder;
	public Configuration config;
	
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
		this.config = new Configuration();
		System.out.println("Current configuration:\n" + config.toString());
	}

	public CryDetectorClassifier(Configuration config) {
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
		ClusteringManager clustering = null;
		System.out.println("#CLUSTERING - START#");
		double[] thresholds = Utils.initializeVector(featureMatrix[0].length, 1);
		for (int tryn = 1; tryn <= maxTriesToFindValidIslandsAndClusters; tryn++) {
			clustering = new ClusteringManager(config, toneUnitCleanedFile);
			clustering.clusterFeatures(featureMatrix, toneUnitCleanedFile.getParentFile(), 1, nClassesDetection, thresholds);

			if (!clustering.highriskclusterfound && tryn == maxTriesToFindValidIslandsAndClusters) {
				System.out.println("A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT!");
				return(null);
			} else if (!clustering.highriskclusterfound) {
				System.out.println(
						"A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT - RETRYING WITH THRESHOLDS REDUCED - RETRY N. "
								+ tryn);
				thresholds = Utils.reduceVectorValues(thresholds, reductionFactor4Retry);
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
		double[][] mod_spec_features_nodelta = Utils.subsetRows(mod_spec_features, numberOfMSFeatures4Detection,
				(numberOfMSFeatures4Detection * 2 - 1)); // mod spec without deltas
		Configuration tempConfig = config.clone();
		System.out.println("\tTransposing matrix for clustering");
		//N VECTORS x N FEATURES 
		double[][] mod_spec_features_nodelta_tr = Utils.traspose(mod_spec_features_nodelta);
		double[] q3 = Utils.columnQ3(mod_spec_features_nodelta_tr);
		System.out.println("\tEnd matrix transposing");

		tempConfig.standardiseFeatures = false;
		ClusteringManager clusteringMS = null;
		for (int tryn = 1; tryn <= maxTriesToFindValidIslandsAndClusters; tryn++) {
			clusteringMS = new ClusteringManager(tempConfig, islandAudio);
			clusteringMS.clusterFeatures(mod_spec_features_nodelta_tr, islandAudio.getParentFile(), nClassesDetection, nClassesDetection, q3);
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
				q3 = Utils.reduceVectorValues(q3, reductionFactor4Retry);
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
		List<double[]> crySegments = SignalConverter.getAnnotatedSignalSegments(islandAudio, clusteringMS.times,clusteringMS.labels); 
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
		short[] reducedsignalMS = SignalConverter.extractAnnotatedSignal(inputAudioReference, times,labels);
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

	public File LSTMDetection(int optimalCluster, File outputFolder) throws Exception {

		System.out.println("Training the LSTM...");
		File modelFolder = dlo.trainLSTM(clusterer.clusters, optimalCluster, config.nhidden, config.nClasses,
				config.minibatch, config.nEpochs);

		System.out.println("Annotating with LSTM...");

		int classifications[] = dlo.annotate(featureMatrix, modelFolder, config.nhidden, config.nClasses,
				config.minibatch, config.nEpochs);
		int nonzero = 0;
		for (int c : classifications) {
			if (c == 1)
				nonzero++;
		}

		File LSTMClusteringFileLab = new File(outputFolder, "LSTM.lab");
		File LSTMClusteringFileLabOutput = new File(outputFolder, outputFolder.getName() + "_LSTM_CTC.lab");

		CorpusCleaner.annotationVector2Lab(classifications, LSTMClusteringFileLab, config.featurewindowsize,
				config.featurewindowshift, 0);
		CorpusCleaner.uniformLabelSegments(LSTMClusteringFileLab, LSTMClusteringFileLabOutput);

		System.out.println("Deleting temp folder " + modelFolder.getName());
		FileUtils.deleteDirectory(modelFolder);
		System.out.println("Training accuracy : " + dlo.trainingAccuracy);
		System.out.println("Evaluation : " + dlo.evaluation);
		System.out.println("Found " + nonzero + " segments");

		return LSTMClusteringFileLabOutput;
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
		double [][] msOverallMatrix = Utils.traspose(ms.modulationSpectrogram);
		ClassificationManager classifierMS = new ClassificationManager(config,unifiedAudioFile);
		File labFile = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", ".lab"));
		//classification clustering 
		classifierMS.clusterFeatures(msOverallMatrix, mainAudioFolder, nClassesClassification, nClassesClassification, entropyEnergyThr, lowestEntropyEnergyThr);
		
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
		msOverallMatrix = Utils.traspose(ms.modulationSpectrogram);
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
		classifierMS.clusterFeatures(msOverallMatrix, mainAudioFolder, nClassesClassification, nClassesClassification, entropyEnergyThr, lowestEntropyEnergyThr);
		
		int signalLengthMS = totalSamples / ModulationSpectrogram.reductionFactor;
		int samplingFreqMS = (int) (af.getSampleRate() / ModulationSpectrogram.reductionFactor);
		double windowShiftMS = ModulationSpectrogram.windowShift;
		
		classifierMS.toLabCTC(labFile, samplingFreqMS, signalLengthMS, windowShiftMS, minimumMSContinuosWindowClassification);
		
		File mod_spec_audio = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_anomalouscry.wav"));
		saveNonEmptyAnnotatedSignal(mod_spec_audio, unifiedAudioFile, classifierMS.times, classifierMS.labels);
		
		
		System.out.println("#CLASSIFICATION FINISHED");
	}

	public static void main(String[] args) throws Exception {

		Configuration config = new Configuration();

		config.maxSilence = 0.5f;
		config.minimumAudioLength = 5f;
		config.energyWindow4Analysis = 0.1f;
		config.pitchWindow4Analysis = 0.1f;
		config.featurewindowsize = 0.3f;
		config.featurewindowshift = 0.1f;
		config.minNFeaturesInCluster = 5;
		config.nClasses = 2;
		config.nhidden = 3;
		config.minibatch = 150;
		config.nEpochs = 2;
		config.standardiseFeatures = true;

		//several short cry impulses
		File audio1 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_01.wav"); //KO //SNR 4.9 duration 16s infant cry contained 0.55s found 0s 

		//long cry followed by short cry impulses with various type of noise in the room (no machine noise) 
		File audio2 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_02.wav"); //SNR 4.67 duration 53s infant cry contained 19s found 1.4 KO 0s

		//short cry imoulses followed by long cry impulses with various type of noise in the room (no machine noise)
		File audio3 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_03.wav"); //SNR 14.95 duration 58s infant cry contained 16s found 1.2s KO 0s
				
		// mostly cry
		File audio4 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_04.wav"); //SNR 8.49 duration 92s infant cry contained 30s found 5.2 KO 0s

		//a lot of overlapping machine beeping
		File audio5 = new File("./test_wave_files/Subintensive-hunger-1child-21_02.wav"); //SNR 16.70 duration 50s infant cry contained 4.1s found 1.69 KO 0s
		
		//lots of adult voices talking aloud
		File audio6 = new File("./test_wave_files/Subintensive-hunger-1child-21_04.wav"); //SNR 6.03 duration 53s infant cry contained 5.4s found 3.48 KO 0s 
		
		//lots of adult voices
		File audio7 = new File("./test_wave_files/Subintensive-pathological_brain-1child#1-23_06.wav"); //SNR 6.8 duration 66s infant cry contained 8s found 5.1s KO 0.5s 

		//poor adult voices overlap telephone and beeps
		File audio8 = new File("./test_wave_files/Subintensive-pathological_brain-2children#1#2-23_05.wav");//SNR 13.47 327s infant cry contained 28s found 1.75s KO 0s

		// cry, noise, and adult voices
		File audio9 = new File("./test_wave_files/Subintensive-unknown-1child-21_01.wav"); //SNR 19.3 duration 152s infant cry contained 4s found 6s KO 2s 

		//no cry present
		File audio10 = new File("./test_wave_files/Subintensive-wakeup-1child-21_05.wav"); //KO //SNR 12.3 duration 56s infant cry contained 0s found 0
		
		//total duration = 1250s = 20.8 min = 1 250 000 energy windows = 100 000 windows of ms (250 ms) extracted every 12ms 
		
		File allAudioToAnalyse [] = {
				audio7,
				//audio8,
				audio1,audio2,audio3,audio4,audio5,audio6,audio9,audio10
				//audio1,audio2,audio3,audio4,audio5,audio6,audio7,audio8,audio9,audio10
		};
		
		/*
		for (File audio:allAudioToAnalyse) {
			System.out.println("\n#######PROCESSING AUDIO FILE :"+audio.getName());
			CryDetectorRefactored cryd = new CryDetectorRefactored(config);
			File output = cryd.detect(audio);
			System.out.println("\n");
		}
		*/
		CryDetectorClassifier cryd = new CryDetectorClassifier(config);
		cryd.classify(allAudioToAnalyse);
	}

}
