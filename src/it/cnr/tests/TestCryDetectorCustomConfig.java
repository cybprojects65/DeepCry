package it.cnr.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.workflow.Configuration;
import it.cnr.workflow.CryDetector;

public class TestCryDetectorCustomConfig {

	public static void main(String[] args) throws Exception{
		
		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1.wav");
		Configuration config = new Configuration();
		
		config.featurewindowsize=0.3f;
		config.featurewindowshift=0.1f;
		
		config.minNFeaturesInCluster = 5;
		
		config.nhidden=3;
		config.minibatch=150;
		config.nEpochs=2;
		
		CryDetector cryd = new CryDetector(config);
		List<File> outputAnnotations = cryd.detect(audio);
		cryd.eval(audio);
		
	}
	
}
