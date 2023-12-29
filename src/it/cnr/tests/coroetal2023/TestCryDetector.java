package it.cnr.tests.coroetal2023;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.workflow.configuration.Configuration;
import it.cnr.workflow.coroetal2023.CryDetectorCoroetAl2023;

public class TestCryDetector {

	public static void main(String[] args) throws Exception{
		
		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1.wav");
		
		
		CryDetectorCoroetAl2023 cryd = new CryDetectorCoroetAl2023();
		List<File> outputAnnotations = cryd.detect(audio);
		cryd.eval(audio);
		
	}
	
}
