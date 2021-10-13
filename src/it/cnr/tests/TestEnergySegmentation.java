package it.cnr.tests;

import java.io.File;

import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;

public class TestEnergySegmentation {

	public static void main(String[] args) throws Exception{

		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1.wav");
		float maxSilence = 0.500f; // s
		float minimumAudioLength = 5; // s
		float energyWindow4Analysis = 0.1f; //s
		float pitchWindow4Analysis = 0.1f; //s
		float featurewindowsize = 1; //s 
		int minNFeaturesInCluster = 2;
		
		FeatureExtractor extractor = new FeatureExtractor();
		
		//works with 44100 kHz audio
		File outputFolder = extractor.separateFilesBasedOnEnergy(audio, maxSilence);
		CorpusCleaner.deleteShortAudio(outputFolder, minimumAudioLength);
		System.out.println("General SNR "+extractor.SNR);
		
		
		File [] outputFolders = extractor.extractEnergyPitchTimeSeriesForAllWaves(outputFolder, energyWindow4Analysis, pitchWindow4Analysis);
		
		for (File fold:outputFolders) {
			
			File [] csvFiles = fold.listFiles();
			for (File csv:csvFiles) {
				if (csv.getName().endsWith(".csv"))
				{
					if (csv.getName().contains("pitch"))
						CorpusCleaner.rawFormat2Lab(csv, new File(csv.getAbsolutePath().replace(".csv", ".lab")), pitchWindow4Analysis);
					else
						CorpusCleaner.rawFormat2Lab(csv, new File(csv.getAbsolutePath().replace(".csv", ".lab")), energyWindow4Analysis);
					
				}
				
			}
			
			//double[][] featureMatrix = extractor.chunkizeTimeSeries(fold, energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize);
			//FeatureClusterer clusterer = new FeatureClusterer();
			//clusterer.clusterFeatures(featureMatrix, fold, minNFeaturesInCluster);
			
			//add clustering to lab conversion
			
			//use another XMeans implementation
			
			//System.exit(0);
		}
		
		
	}

}
