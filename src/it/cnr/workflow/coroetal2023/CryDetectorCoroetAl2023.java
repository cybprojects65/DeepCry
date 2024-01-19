package it.cnr.workflow.coroetal2023;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;

import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.evaluation.Evaluator;
import it.cnr.evaluation.Range;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.utilities.SignalProcessing;

public class CryDetectorCoroetAl2023 {

	public MultiKMeans clusterer = new MultiKMeans();
	public DeepLearningManager dlo = new DeepLearningManager();
	public double[][] featureMatrix = null;
	public List<File> extractedWaveFiles = null;
	public List<File> annotationFiles = null;
	public double SNR = 0;
	
	public double maxCryAnnotationPerc = 0.98;

	public CryDetectorCoroetAl2023() {
		this.config = new WorkflowConfiguration();
		System.out.println("Current configuration:\n" + config.toString());
	}

	public WorkflowConfiguration config;

	public CryDetectorCoroetAl2023(WorkflowConfiguration config) {
		this.config = config.clone();
		System.out.println("Current configuration:\n" + config.toString());
	}

	public List<File> detect(File audio44100Hz) throws Exception {

		// step 1 - separate input file into subfiles based on energy
		File[] subfilefeaturefolders = segmentInputFile(audio44100Hz);

		annotationFiles = new ArrayList<>();

		for (File featuresfolder : subfilefeaturefolders) {
			// step 2 - run cluster analysis on each subfolder and detect optimal cluster
			System.out.println("Step 2 - Cluster analysis and optimal cluster detection");
			int optimalCluster = clusterAnalysis(featuresfolder);
			if (optimalCluster > -1) {
				System.out.println("Optimal cluster size for " + featuresfolder.getName() + " is "
						+ clusterer.clusters.get(optimalCluster).size());
				// step 3 - train an LSTM and apply it to recognize further vectors
				System.out.println("Step 3 - Training and applying LSTM");
				File annotation = LSTMDetection(optimalCluster, featuresfolder);
				annotationFiles.add(annotation);
				System.out.println("Step F - Wokflow complete for file " + featuresfolder.getName());
			}else {
				System.out.println("----------------------WARNING----------------------------");
				System.out.println("Segment "+featuresfolder.getName()+" MAY not contain CRY!");
				System.out.println("---------------------------------------------------------");
				Thread.sleep(3000);
				File LSTMClusteringFileLabOutput = new File(featuresfolder, featuresfolder.getName() + "_LSTM_CTC.lab");
				LSTMClusteringFileLabOutput.createNewFile();
				annotationFiles.add(LSTMClusteringFileLabOutput);
			}
		}

		System.out.println("Cry detection finished for file " + audio44100Hz.getName());
		
		//check annotation consistency
		//reviseAnnotations();
		
		return annotationFiles;

	}

	public void reviseAnnotations() throws Exception {
		
		
		for (File annotation:annotationFiles) {
			Evaluator eval = new Evaluator();
			String waveFile = annotation.getName();
			waveFile = waveFile.substring(0,waveFile.indexOf(".wav"))+".wav";
			File correspondentWaveFile = new File(annotation.getParentFile().getParentFile(),waveFile);
			AudioBits ab = new AudioBits(correspondentWaveFile);
			float sfrequency = ab.getAudioFormat().getSampleRate();
			int nsamples = ab.getShortVectorAudio().length;
			ab.ais.close();
			double durationinsec = SignalProcessing.samplesToTime(nsamples, (double)sfrequency);
			
			Range [] ranges = eval.annotationToRanges(annotation);
			double annotatedDuration = 0;
			for (Range r:ranges) {
				double d = r.t1-r.t0;
				annotatedDuration+=d;
			}
			
			double percDur = annotatedDuration/durationinsec;
			System.out.println("Perc. annotation of "+annotation.getName()+"="+percDur+" ("+annotatedDuration+" vs "+durationinsec+")");
			if (percDur > maxCryAnnotationPerc) {
				
				System.out.println("----------------------WARNING----------------------------");
				System.out.println("Emptying annotation for segment "+annotation.getName()+" because it overestimated cry segments!");
				System.out.println("---------------------------------------------------------");
				Thread.sleep(3000);
				
				FileUtils.forceDelete(annotation);
				annotation.createNewFile();
			}
			
			
		}//end for
		
	}
	
	
	public File[] segmentInputFile(File audio44100Hz) throws Exception {

		FeatureExtractor extractor = new FeatureExtractor();

		// works with 44100 kHz audio
		File outputFolder = extractor.separateFilesBasedOnEnergy(audio44100Hz, config.maxSilence);
		SNR = extractor.SNR;
		CorpusCleaner.deleteShortAudio(outputFolder, config.minimumAudioLength);
		System.out.println("General SNR " + extractor.SNR);

		File[] outputFolders = extractor.extractEnergyPitchTimeSeriesForAllWaves(outputFolder,
				config.energyWindow4Analysis, config.pitchWindow4Analysis);

		extractedWaveFiles = new ArrayList<File>();

		for (File fold : outputFolders) {
			String originalWaveFile = fold.getAbsolutePath();
			originalWaveFile = originalWaveFile.substring(0, originalWaveFile.indexOf(".wav")) + ".wav";
			extractedWaveFiles.add(new File(originalWaveFile));

			File[] csvFiles = fold.listFiles();
			for (File csv : csvFiles) {
				if (csv.getName().endsWith(".csv")) {
					if (csv.getName().contains("pitch"))
						CorpusCleaner.rawFormat2Lab(csv, new File(csv.getAbsolutePath().replace(".csv", ".lab")),
								config.pitchWindow4Analysis);
					else
						CorpusCleaner.rawFormat2Lab(csv, new File(csv.getAbsolutePath().replace(".csv", ".lab")),
								config.energyWindow4Analysis);

				}

			}

		}

		return outputFolders;

	}

	public int clusterAnalysis(File featureFolder) throws Exception {
		clusterer = new MultiKMeans();
		FeatureExtractor extractor = new FeatureExtractor();
		// retrieve the matrix of the features
		featureMatrix = extractor.chunkizeTimeSeries(featureFolder, config.energyWindow4Analysis,
				config.pitchWindow4Analysis, config.featurewindowsize, config.featurewindowshift);
		File clusteringFile = clusterer.clusterFeatures(featureMatrix, featureFolder, config.minNFeaturesInCluster);
		// save the model
		File outputClusteringModelFile = new File(featureFolder, "model_clustering.bin");
		clusterer.save(outputClusteringModelFile);
		// save the clustering file
		File clusteringFileLab = new File(clusteringFile.getAbsolutePath().replace(".csv", ".lab"));
		CorpusCleaner.clusteringFormat2Lab(clusteringFile, clusteringFileLab, config.featurewindowsize,
				config.featurewindowshift);
		dlo = new DeepLearningManager();
		// detect the most cry-like cluster and the shape of the time series
		//NOTE: SET THIS TO HMM2 TO DISCARD NON CRY SEGMENTS
		int optimalCluster = dlo.calcMaxLikelyhoodClusterHMM2(clusterer.clusters, config.energyWindow4Analysis,
				config.pitchWindow4Analysis, config.featurewindowsize, config.featurewindowshift,
				config.minNFeaturesInCluster, clusterer.Kstar);

		return optimalCluster;

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

	public double accuracy;
	public double precision;
	public double recall;
	public double f1;
	Evaluator evaluator;

	public void eval(File audio) throws Exception {
		List<File> waveFiles = extractedWaveFiles;
		List<File> goldenAnnotations = new ArrayList<>();
		for (File w : waveFiles) {
			goldenAnnotations.add(new File(audio.getParentFile(), w.getName().replace(".wav", ".lab")));

		}

		evaluator = new Evaluator();
		evaluator.evaluate(annotationFiles, goldenAnnotations, waveFiles);
		accuracy = evaluator.accuracy;
		precision = evaluator.precision;
		recall = evaluator.recall;
		f1 = evaluator.f1;

	}
}
