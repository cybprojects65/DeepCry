package it.cnr.hmm;

import java.io.File;

import it.cnr.features.CorpusCleaner;
import it.cnr.features.EnergyPitchFeatureExtractor;
import it.cnr.workflow.configuration.WorkflowConfiguration;

public class CreateReferenceFeatures {

	public static void main(String[] args) throws Exception{

		EnergyPitchFeatureExtractor extractor = new EnergyPitchFeatureExtractor();
		File outputFolder = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\Cry only");
		WorkflowConfiguration config = new WorkflowConfiguration();
		
		File [] outputFolders = extractor.extractEnergyPitchTimeSeriesForAllWaves(outputFolder, config.energyWindow4Analysis, config.pitchWindow4Analysis);
		
		
		
		
	}
	
}
