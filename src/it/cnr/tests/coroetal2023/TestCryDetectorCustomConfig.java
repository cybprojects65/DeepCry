package it.cnr.tests.coroetal2023;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.workflow.configuration.Configuration;
import it.cnr.workflow.coroetal2023.CryDetectorCoroetAl2023;

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
		
		CryDetectorCoroetAl2023 cryd = new CryDetectorCoroetAl2023(config);
		List<File> outputAnnotations = cryd.detect(audio);
		cryd.eval(audio);
		
	}
	
}
