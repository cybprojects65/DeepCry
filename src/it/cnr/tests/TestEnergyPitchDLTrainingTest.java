package it.cnr.tests;

import java.io.File;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;

import it.cnr.clustering.MultiKMeans;
import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;

public class TestEnergyPitchDLTrainingTest {

	public static void main(String[] args) throws Exception{

		File fold = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\");
		float energyWindow4Analysis = 0.1f; //s
		float pitchWindow4Analysis = 0.1f; //s
		float featurewindowsize = 0.4f; //s 
		float featurewindowshift = 0.1f; //s 
		int minNFeaturesInCluster = 3;
		int nClasses = 2;
		
		int nhidden = 3;//3;
		int minibatch = 50;  //90
		int nEpochs = 10;
		
		/*
		int nhidden = 25;
		int minibatch = 90;
		int nEpochs = 10;
		*/
		//TODO: 
		//0 - add reference to training set - KO
		//1 - produce complete workflow
		//2 - annotate audio 3 too
		//3 - evaluate performance wrt my annotations
		//3 - reduce feature windowsize to 0.3 OK - worse
		
		File outputClusteringFile = new File("temp_clustering.bin");
		MultiKMeans clusterer = MultiKMeans.load(outputClusteringFile);
		
		FeatureExtractor extractor = new FeatureExtractor();
		double[][] featureMatrix = extractor.chunkizeTimeSeries(fold, energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize,featurewindowshift);
		
		DeepLearningManager dlo = new DeepLearningManager();
		int optimalCluster = dlo.calcMaxLikelyhoodClusterHMM(clusterer.clusters, energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift, minNFeaturesInCluster, clusterer.Kstar);
		
		HashMap<Integer, List<double[]>> references = dlo.getReferenceVectors(energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift);
		
		List<double[]> optcluster = clusterer.clusters.get(optimalCluster);
		
		for (Integer r:references.keySet()) {
			
			//optcluster.addAll(references.get(r));
			
		}
		
		System.out.println("Optimal cluster size is "+clusterer.clusters.get(optimalCluster).size());
		
		//dlo.buildTrainingSet(clusterer.clusters,optimalCluster);
		File modelFolder = dlo.trainLSTM(clusterer.clusters, optimalCluster,
				nhidden,
				nClasses,
				minibatch,
				nEpochs);
		
		
		
		int classifications[] = dlo.annotate(featureMatrix, modelFolder, nhidden,	nClasses,	minibatch,	nEpochs);
		int nonzero = 0;
		for (int c:classifications) {
			if (c==1)
				nonzero++;
		}
		
		File LSTMClusteringFileLab = new File(fold, "LSTM.lab");
		
		CorpusCleaner.annotationVector2Lab(classifications,LSTMClusteringFileLab, featurewindowsize,featurewindowshift,0);
		
		System.out.println("Deleting temp folder "+modelFolder.getName());
		FileUtils.deleteDirectory(modelFolder);
		//OK - train an LSTM
		System.out.println("Training accuracy : "+dlo.trainingAccuracy);
		System.out.println("Evaluation : "+dlo.evaluation);
		System.out.println("Found "+nonzero+" segments");
	}
		
		

}
