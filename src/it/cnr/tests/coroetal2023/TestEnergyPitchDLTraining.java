package it.cnr.tests.coroetal2023;

import java.io.File;
import java.util.HashMap;
import java.util.List;

import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.FeatureExtractor;

public class TestEnergyPitchDLTraining {

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
		
		//dlo.buildTrainingSet(clusterer.clusters,optimalCluster);
		dlo.trainLSTM(clusterer.clusters, optimalCluster,
				nhidden,
				nClasses,
				minibatch,
				nEpochs);
		
		System.out.println("Training accuracy : "+dlo.trainingAccuracy);
		
		//OK - train an LSTM
		
	}
		
		

}
