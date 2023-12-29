package it.cnr.hmm;

import java.io.File;

import it.cnr.features.CorpusCleaner;
import it.cnr.features.FeatureExtractor;
import it.cnr.workflow.configuration.Configuration;

public class CreateReferenceFeatures {

	public static void main(String[] args) throws Exception{

		FeatureExtractor extractor = new FeatureExtractor();
		File outputFolder = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\Cry only");
		Configuration config = new Configuration();
		
		File [] outputFolders = extractor.extractEnergyPitchTimeSeriesForAllWaves(outputFolder, config.energyWindow4Analysis, config.pitchWindow4Analysis);
		
		
		
		
	}
	
}
