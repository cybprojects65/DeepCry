package it.cnr.tests.coroetal2024;

import java.io.File;

import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.early.CryDetectorClassifier;
import it.cnr.workflow.coroetal2024.staging.AnomalousCryDetector;

public class TestAnomalousCryDetector {

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
		config.nhidden = 3;
		config.minibatch = 150;
		config.nEpochs = 2;
		config.standardiseFeatures = true;

		//several short cry impulses
		File audio1 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_01.wav"); //KO //SNR 4.9 duration 16s infant cry contained 0.55s found 0s 

		//long cry followed by short cry impulses with various type of noise in the room (no machine noise) 
		File audio2 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_02.wav"); //SNR 4.67 duration 53s infant cry contained 19s found 1.4 KO 0s

		//short cry imoulses followed by long cry impulses with various type of noise in the room (no machine noise)
		File audio3 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_03.wav"); //SNR 14.95 duration 58s infant cry contained 16s found 1.2s KO 0s
				
		// mostly cry
		File audio4 = new File("./test_wave_files/Subintensive-annoyed-1child#2-23_04.wav"); //SNR 8.49 duration 92s infant cry contained 30s found 5.2 KO 0s

		//a lot of overlapping machine beeping
		File audio5 = new File("./test_wave_files/Subintensive-hunger-1child-21_02.wav"); //SNR 16.70 duration 50s infant cry contained 4.1s found 1.69 KO 0s
		
		//lots of adult voices talking aloud
		File audio6 = new File("./test_wave_files/Subintensive-hunger-1child-21_04.wav"); //SNR 6.03 duration 53s infant cry contained 5.4s found 3.48 KO 0s 
		
		//lots of adult voices
		File audio7 = new File("./test_wave_files/Subintensive-pathological_brain-1child#1-23_06.wav"); //SNR 6.8 duration 66s infant cry contained 8s found 5.1s KO 0.5s 

		//poor adult voices overlap telephone and beeps
		File audio8 = new File("./test_wave_files/Subintensive-pathological_brain-2children#1#2-23_05.wav");//SNR 13.47 327s infant cry contained 28s found 1.75s KO 0s

		// cry, noise, and adult voices
		File audio9 = new File("./test_wave_files/Subintensive-unknown-1child-21_01.wav"); //SNR 19.3 duration 152s infant cry contained 4s found 6s KO 2s 

		//no cry present
		File audio10 = new File("./test_wave_files/Subintensive-wakeup-1child-21_05.wav"); //KO //SNR 12.3 duration 56s infant cry contained 0s found 0
		
		//total duration = 1250s = 20.8 min = 1 250 000 energy windows = 100 000 windows of ms (250 ms) extracted every 12ms 
		
		File allAudioToAnalyse [] = {
				audio7
				,
				audio8,
				audio1,audio2,audio3,audio4,audio5,audio6,
				audio9 
				,audio10
		};
		
		
		AnomalousCryDetector cryd = new AnomalousCryDetector(config);
		cryd.run(allAudioToAnalyse);
	}

}
