package it.cnr.tests.coroetal2023;

import java.io.File;

import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;

public class TestEnergyPitchDLClassify {

	public static void main(String[] args) throws Exception{

		File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\");
		float energyWindow4Analysis = 0.1f; //s
		float pitchWindow4Analysis = 0.1f; //s
		float featurewindowsize = 0.4f; //s 
		float featurewindowshift = 0.1f; //s 
		int minNFeaturesInCluster = 3;
		int nhidden = 3;
		int nClasses = 2;
		int minibatch = 50;
		int nEpochs = 10;
		
		File outputClusteringFile = new File("temp_clustering.bin");
		MultiKMeans clusterer = MultiKMeans.load(outputClusteringFile);
		
		FeatureExtractor extractor = new FeatureExtractor();
		double[][] featureMatrix = extractor.chunkizeTimeSeries(fold, energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize,featurewindowshift);
		DeepLearningManager dlo = new DeepLearningManager();
		int optimalCluster = dlo.calcMaxLikelyhoodClusterHMM(clusterer.clusters, energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift, minNFeaturesInCluster, clusterer.Kstar);
		
		File modelFolder = new File("tempb6889ec4-55cc-4fa6-a35f-1563fda89558");
		int classifications[] = dlo.annotate(featureMatrix, modelFolder, nhidden,	nClasses,	minibatch,	nEpochs);
		
		File LSTMClusteringFileLab = new File(fold, "LSTM.lab");
		
		CorpusCleaner.annotationVector2Lab(classifications,LSTMClusteringFileLab, featurewindowsize,featurewindowshift,0);
		
		
		//OK - annotate file through the LSTM
		//ADDITIONAL: ADD spectral features
		
	}
		
		

}
