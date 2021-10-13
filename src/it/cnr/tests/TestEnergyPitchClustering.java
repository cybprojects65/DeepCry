package it.cnr.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;
import it.cnr.workflow.Configuration;

public class TestEnergyPitchClustering {

	public static void main(String[] args) throws Exception{

		//File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\");
		//File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB1\\energysegmentation_0.5\\audio_segment_18.wav_energysegmentation_e0.1_p0.1");
		//File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB1\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1");
		File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB1\\energysegmentation_0.5\\audio_segment_1.wav_energysegmentation_e0.1_p0.1");
		//File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB1\\energysegmentation_0.5\\audio_segment_14.wav_energysegmentation_e0.1_p0.1");
		Configuration c = new Configuration(); 
		
		MultiKMeans clusterer = new MultiKMeans();
		FeatureExtractor extractor = new FeatureExtractor();
		
		double[][] featureMatrix = extractor.chunkizeTimeSeries(fold, c.energyWindow4Analysis, c.pitchWindow4Analysis, c.featurewindowsize,c.featurewindowshift);
		File clusteringFile = clusterer.clusterFeatures(featureMatrix, fold, c.minNFeaturesInCluster);
		
		DeepLearningManager dlo = new DeepLearningManager();
		int optimalCluster = dlo.calcMaxLikelyhoodClusterHMM2(clusterer.clusters, c.energyWindow4Analysis, c.pitchWindow4Analysis, c.featurewindowsize, c.featurewindowshift, c.minNFeaturesInCluster, clusterer.Kstar);
		
		File clusteringFileLab = new File(clusteringFile.getAbsolutePath().replace(".csv", ".lab"));
		CorpusCleaner.clusteringFormat2Lab(clusteringFile, clusteringFileLab, c.featurewindowsize,c.featurewindowshift);
		List<Integer> goodclustersList = new ArrayList<>(); 
		goodclustersList.add(optimalCluster);
		File goodClusteringFileLab = new File(clusteringFile.getAbsolutePath().replace(".csv", "_goodclusters.lab"));
		CorpusCleaner.clusteringFormat2Lab(clusteringFile, goodClusteringFileLab, c.featurewindowsize,c.featurewindowshift, goodclustersList);
		
	}
		
		

}
