package it.cnr.tests.coroetal2024;

import java.io.File;

import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.CryDetectorClassifier;

public class TestModulationSpectrogram {

	public static void main(String[] args) throws Exception {

		boolean saturateMagnitudeDBs = true;
		boolean addDeltas = false;
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
		CryDetectorClassifier cdc = new CryDetectorClassifier(config);
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-pathological_brain-1child#1-23_06_processing\\Subintensive-pathological_brain-1child#1-23_06_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-pathological_brain-2children#1#2-23_05_processing\\Subintensive-pathological_brain-2children#1#2-23_05_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-annoyed-1child#2-23_04_processing\\Subintensive-annoyed-1child#2-23_04_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-annoyed-1child#2-23_03_processing\\Subintensive-annoyed-1child#2-23_03_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-annoyed-1child#2-23_02_processing\\Subintensive-annoyed-1child#2-23_02_shortened_islands.wav");
		
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-wakeup-1child-21_05_processing\\Subintensive-wakeup-1child-21_05_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-unknown-1child-21_01_processing\\Subintensive-unknown-1child-21_01_shortened_islands.wav");
		//File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-annoyed-1child#2-23_01_processing\\Subintensive-annoyed-1child#2-23_01_shortened.wav");
		
		File testAudio = new File("C:\\Users\\Utente\\eclipse-workspace-DL\\DeepCry\\test_wave_files\\Subintensive-pathological_brain-2children#1#2-23_05_processing\\Subintensive-pathological_brain-2children#1#2-23_05_shortened.wav");
		
		
		File ms_output = new File(testAudio.getAbsolutePath().replace(".wav", "_mod_spect_forTesting.csv"));
		
		
		
		System.out.println("#MODULATION SPECTROGRAM - START FOR FILE "+testAudio.getAbsolutePath());
		ModulationSpectrogram ms = new ModulationSpectrogram();
		ms.calcMS(testAudio, ms_output, saturateMagnitudeDBs, addDeltas, cdc.numberOfMSFeatures4Classification, cdc.maxFrequencyForClassification);
		System.out.println("#MODULATION SPECTROGRAM - END FOR FILE "+testAudio.getAbsolutePath());

	}

}
