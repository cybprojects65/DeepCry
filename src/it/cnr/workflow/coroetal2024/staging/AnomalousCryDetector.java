package it.cnr.workflow.coroetal2024.staging;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.sound.sampled.AudioFormat;

import it.cnr.anomaly.DLAnomalyDetection;
import it.cnr.anomaly.MSEnergyEntropyAnomalyDetection;
import it.cnr.clustering.ClassificationManager;
import it.cnr.clustering.DetectionManager;
import it.cnr.clustering.EnergyPitchFilterManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.EnergyPitchFeatureExtractor;
import it.cnr.features.FeatureSet;
import it.cnr.features.IslandDetector;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.AnomalousCryConfiguration;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsObjects;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class AnomalousCryDetector {

	int signalLengthUnifiedFileInMS; 
    int samplingFreqUnifiedFile;
    double[] classificationTimeRanges;
    String[] classificationLabels;
    
	
	
	public AnomalousCryDetector() {
		
	}

	public File run(File allAudiotoAnalyse[]) throws Exception {

		if(!AnomalousCryConfiguration.skippreprocessing) {
		System.out.println("#################################################");
		System.out.println("######### PROCESSING INDIVIDUAL FILES ###########");
		System.out.println("#################################################");
		for (File audio : allAudiotoAnalyse) {
			System.out.println("\n#######PROCESSING AUDIO FILE :" + audio.getName());
			// 1 - tone unit segmentation
			System.out.println("#1- TONE UNIT SEGMENTATION - START#");
			File toneUnitCleanedFile = toneUnitSegmentation(audio);
			System.out.println("#1- TONE UNIT SEGMENTATION - END#");
			// 2 - energy/pitch feature extration at 100ms
			System.out.println("#2 - ENERGY/PITCH FEATURE EXTRACTION - START#");
			EnergyPitchFeatureExtractor featurer = new EnergyPitchFeatureExtractor(AnomalousCryConfiguration.window4EnergyPitchFeatures);
			// one row for each frame
			double featureMatrix[][] = featurer.extractFeatureMatrix(toneUnitCleanedFile, true);
			System.out.println("#2 - ENERGY/PITCH FEATURE EXTRACTION - END#");
			// 3 - energy/pitch clustering
			System.out.println("#3 - HIGH ENERGY/PITCH SEGMENT DETECTION - START#");
			DetectionManager clustering = clusterFeaturesForDetection(featureMatrix, toneUnitCleanedFile);
			System.out.println("#3 - HIGH ENERGY/PITCH SEGMENT DETECTION - END#");
			// 4 - energy island identification
			if (clustering != null) {
				System.out.println("#4 - ENERGY ISLAND DETECTION - START#");
				IslandDetector islandDetector = new IslandDetector(clustering, AnomalousCryConfiguration.window4EnergyPitchFeatures, toneUnitCleanedFile);
				islandDetector.detectIslands();
				System.out.println("#4 - ENERGY ISLAND DETECTION - END#");
			}
			
		}
		}
		
		System.out.println("#################################################");
		System.out.println("######### PROCESSING MERGED FILE ################");
		System.out.println("#################################################");
		
		System.out.println("#5 - SIGNAL ISLAND MERGING - START#");
		File unifiedAudioFile = mergeIslands(allAudiotoAnalyse);
		System.out.println("#5 - SIGNAL ISLAND MERGING - START#");
		File mscacheFile = new File(unifiedAudioFile.getParentFile(),"ms.bin");
		File msCSVFile = new File(unifiedAudioFile.getParentFile(),"ms.csv");
		
		System.out.println("#6 - CALCULATING LOW FREQUENCY MODULATION SPECTROGRAM - START#");
		double [][] msOverallMatrix = null;
		if (!AnomalousCryConfiguration.skippreprocessing) {
			System.out.println("Unified file is "+unifiedAudioFile.getAbsolutePath());
			ModulationSpectrogram ms = new ModulationSpectrogram();
			boolean saturateMagnitudeDBs = true;
			boolean addDeltas = false;
			ms.calcMS(unifiedAudioFile, msCSVFile, saturateMagnitudeDBs, addDeltas,AnomalousCryConfiguration.numberOfMSFeatures,AnomalousCryConfiguration.maxMSFrequency);
			System.out.println("#6 - CALCULATING LOW FREQUENCY MODULATION SPECTROGRAM - END#");
			//FROM rows = features ; columns = observations -> rows = observations ; columns = features
			msOverallMatrix = UtilsVectorMatrix.traspose(ms.modulationSpectrogram);
			UtilsObjects.saveObject(mscacheFile, msOverallMatrix);
		}else {
			msOverallMatrix = (double [][]) UtilsObjects.loadObject(mscacheFile);
		}
		
		System.out.println("#7 - CLASSIFICATION OF LOW FREQUENCY MODULATION SPECTROGRAM - START#");
		ClassificationManager classifierMS = null;
		if (!AnomalousCryConfiguration.useDLClassification) {
			MSEnergyEntropyAnomalyDetection classifierMSEE = new MSEnergyEntropyAnomalyDetection(unifiedAudioFile);
			classifierMSEE.classifyFeatures(msOverallMatrix,AnomalousCryConfiguration.lowEnergy, AnomalousCryConfiguration.highEnergy, AnomalousCryConfiguration.lowEntropy, AnomalousCryConfiguration.highEntropy);
			classifierMS = classifierMSEE; 
		}else {
			DLAnomalyDetection classifierMSDL = new DLAnomalyDetection(unifiedAudioFile);
			classifierMSDL.classifyFeatures(msOverallMatrix);
			classifierMS = classifierMSDL;
		}
		System.out.println("#7 - CLASSIFICATION OF LOW FREQUENCY MODULATION SPECTROGRAM - END#");
		
		System.out.println("#8 - ANNOTATE SIGNAL WITH THE RAW CLASSIFIED FEATURES - START#");
		classifierMS = annotateSignalwithMSClassification( unifiedAudioFile,  classifierMS);
		System.out.println("#8 - ANNOTATE SIGNAL WITH THE RAW CLASSIFIED FEATURES - END#");
		
		if (classifierMS.highriskclusterfound) {
			System.out.println("#9 - EXTEND THE CLASSIFIED SEGMENTS - START#");
			extendAnnotationThroughEnergyIslands(msOverallMatrix, classifierMS.times, classifierMS.labels, unifiedAudioFile);
			System.out.println("#9 - EXTEND THE CLASSIFIED SEGMENTS - END#");
			System.out.println("#10 - SAVING DETECTED ANOMALOUS CRY IN A SEPARATE FILE - START#");
			File mod_spec_audio = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_anomalouscry.wav"));
			saveNonEmptyAnnotatedSignal(mod_spec_audio, unifiedAudioFile, classificationTimeRanges, classificationLabels);
			System.out.println("#10 - SAVING DETECTED ANOMALOUS CRY IN A SEPARATE FILE - END#");
		}else {
			System.out.println("#####EXIT : NO ANOMALY DETECTED#");
			System.err.println("#####EXIT : NO ANOMALY DETECTED#");
		}
		
		System.out.println("#PROCESS END.#");
		return unifiedAudioFile;
		
	}

	public ClassificationManager annotateSignalwithMSClassification(File unifiedAudioFile, ClassificationManager classifierMS) throws Exception {
		
		File rawLabFile = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", "_msclassification.lab"));
		double windowShiftMS = ModulationSpectrogram.windowShift;
		AudioBits ab = new AudioBits(unifiedAudioFile);
		short [] signal = ab.getShortVectorAudio();
		int totalSamples = signal.length;
		ab.ais.close();
		int signalLengthMS = totalSamples / ModulationSpectrogram.reductionFactor;
		int samplingFreqMS = (int) (ab.getAudioFormat().getSampleRate() / ModulationSpectrogram.reductionFactor);
		classifierMS.toLabCTC(rawLabFile, samplingFreqMS, signalLengthMS, windowShiftMS,AnomalousCryConfiguration.minimumContinuousAnomalousSegment4Classification);
		return classifierMS;
	}
	
	public static void saveNonEmptyAnnotatedSignal(File outputAudiofile, File inputAudioReference, double [] times, String [] labels) throws Exception{
		short[] reducedsignalMS = SignalProcessing.extractAnnotatedSignal(inputAudioReference, times,labels,AnomalousCryConfiguration.silenceSecondsToAddBetweenAnomalousCrySamples4Reporting);
		System.out.println("Saving annotated audio segments");
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(reducedsignalMS, outputAudiofile,
			new AudioBits(inputAudioReference).getAudioFormat());
	}
	
	public File extendAnnotationThroughEnergyIslands(double [][] msOverallMatrix, double classificationtimes[], String [] classificationlabels, File  unifiedAudioFile) throws Exception{
		
		FeatureSet fs = FeatureSet.fromModulationSpectrogram(msOverallMatrix, classificationtimes, classificationlabels, ModulationSpectrogram.windowShift);
		IslandDetector id = new IslandDetector(AnomalousCryConfiguration.window4EnergyPitchFeatures, unifiedAudioFile);
		List<double []> times = id.extendIntervalsThroughEnergyIslands(fs.timeIntervals);
		
		List<Double> timesN = new ArrayList<>();
		List<String> labelsN = new ArrayList<>();
		timesN.add(0d);
		labelsN.add("");
		  
		for (double [] tint : times) {
			  
			  double t0 = tint[0];
			  double t1 = tint[1];
			  String l = "A";
			  timesN.add(t0);
			  labelsN.add(l);
			  timesN.add(t1);
			  labelsN.add("");
			  
		  }
		  
		  classificationTimeRanges = new double[timesN.size()];
		  for (Integer h = 0; h <classificationTimeRanges.length;h++) {
			  classificationTimeRanges[h] = timesN.get(h);
		  }

		  classificationLabels = new String[labelsN.size()];

		  for (Integer h = 0; h <classificationLabels.length;h++) {
			  classificationLabels[h] = labelsN.get(h);
		  }
		  
		  File labeloutput = new File(unifiedAudioFile.getAbsolutePath().replace(".wav", ".lab"));
		  CorpusCleaner.vector2LabFile(classificationLabels, classificationTimeRanges, labeloutput);
		  return labeloutput;
	}
	
	public File mergeIslands(File audioFiles[]) throws Exception{

		AudioFormat af = null;
		int totalSamples = 0;
		List<short[]> signals = new ArrayList<>();
		
		File unifiedAudioFile = null; 
		
		for (File audio : audioFiles) {
			
			File audioFolder = new File(audio.getAbsolutePath().replace(".wav", "_processing"));
			File[] allfolderfiles = audioFolder.listFiles();
			File modspecfile = null;
			for (File f : allfolderfiles) {
				if (f.getName().endsWith("_islands.wav")) {
					modspecfile = f;
					break;
				}
			}
			if (modspecfile == null)
				continue;

			System.out.println("#ADDING AUDIO FILE " + audio.getName());
			AudioBits ab = new AudioBits(modspecfile);
			if (af == null) {
				af = new AudioBits(modspecfile).getAudioFormat();
				unifiedAudioFile = new File (audio.getParentFile(),"islands_merged.wav");
			}
			short[] sig = ab.getShortVectorAudio();
			ab.ais.close();
			totalSamples += sig.length;
			signals.add(sig);
			
			short silence[] = SignalProcessing.silence(AnomalousCryConfiguration.maxSilence4ToneUnitSelection, af.getSampleRate());
			totalSamples += silence.length;
			//add silence
			signals.add(silence);
			
		}

		short[] totalFile = new short[totalSamples];
		int sam = 0;
		for (short[] sig : signals) {
			for (int i = 0; i < sig.length; i++) {
				totalFile[sam] = sig[i];
				sam++;
			}
		}

		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(totalFile, unifiedAudioFile, af);
		System.out.println("Unifying file to " + unifiedAudioFile.getAbsolutePath() + " done");
		
		signalLengthUnifiedFileInMS = totalSamples / ModulationSpectrogram.reductionFactor;
		samplingFreqUnifiedFile = (int) (af.getSampleRate() / ModulationSpectrogram.reductionFactor);
		
		return unifiedAudioFile;
	}

	public DetectionManager clusterFeaturesForDetection(double[][] featureMatrix, File toneUnitCleanedFile)
			throws Exception {
		EnergyPitchFilterManager clustering = null;

		double[] thresholds = UtilsVectorMatrix.initializeVector(featureMatrix[0].length,AnomalousCryConfiguration.minimumScaledEnergyPitch);
		for (int tryn = 1; tryn <= AnomalousCryConfiguration.maxTriesToFindValidIslandsAndClusters; tryn++) {

			clustering = new EnergyPitchFilterManager(toneUnitCleanedFile);
			clustering.detectHighValuedFeatures(featureMatrix, thresholds);

			if (!clustering.highriskclusterfound && tryn == AnomalousCryConfiguration.maxTriesToFindValidIslandsAndClusters) {
				System.out.println("A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT!");
				return (null);
			} else if (!clustering.highriskclusterfound) {
				System.out.println(
						"A HIGH-VALUED CLUSTER WAS NOT FOUND - INFANT CRY POSSIBLY ABSENT - RETRYING WITH THRESHOLDS REDUCED - RETRY N. "
								+ tryn);
				thresholds = UtilsVectorMatrix.reduceVectorValues(thresholds, AnomalousCryConfiguration.reductionFactor4Retry);
				System.out.println("New thresholds: " + Arrays.toString(thresholds));
			} else {
				System.out.println("A HIGH-VALUED CLUSTER WAS FOUND");
				break;
			}

		}

		return clustering;
	}

	public File toneUnitSegmentation(File audio) throws Exception {
		EnergyPitchFeatureExtractor extractor = new EnergyPitchFeatureExtractor(AnomalousCryConfiguration.window4EnergyPitchFeatures);
		// works with 44100 kHz audio
		File outputFolder = extractor.separateFilesBasedOnEnergy(audio, (float) AnomalousCryConfiguration.maxSilence4ToneUnitSelection);
		CorpusCleaner.deleteShortAudio(outputFolder,AnomalousCryConfiguration.minimumAudioLength4ToneUnitSelection);
		System.out.println("SNR " + extractor.SNR);
		File outputCleanedFile = CorpusCleaner.unifyToneUnitsIntoOneFile(outputFolder, audio);
		return outputCleanedFile;
	}

}
