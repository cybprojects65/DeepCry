package it.cnr.tests.coroetal2024;

import java.io.File;

import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.CryDetectorClassifier;

public class TestInfantCryClassifierExtender {

	public static void main(String[] args) throws Exception {

		WorkflowConfiguration config = new WorkflowConfiguration();

		config.maxSilence = 0.5f;
		config.minimumAudioLength = 5f;
		config.energyWindow4Analysis = 0.1f;
		config.pitchWindow4Analysis = 0.1f;
		config.featurewindowsize = 0.3f;
		config.featurewindowshift = 0.1f;
		config.minNFeaturesInCluster = 5;
		config.nClasses = 2;
		config.nhidden = 2; //2
		config.minibatch = 1000;
		config.nEpochs = 1;
		config.standardiseFeatures = true;

		CryDetectorClassifier cryd = new CryDetectorClassifier(config);
		File msmatrix = new File("./test_wave_files/islands_unified.bin");
		File anomalousCryFeatures = new File("./test_wave_files/islands_unified_featureset.bin");
		File unifiedFile = new File("./test_wave_files/islands_unified.wav");
		
		//cryd.extendClassification(msmatrix,anomalousCryFeatures,unifiedFile);
		cryd.extendClassificationIslands(msmatrix,anomalousCryFeatures,unifiedFile);
		
	}

}
