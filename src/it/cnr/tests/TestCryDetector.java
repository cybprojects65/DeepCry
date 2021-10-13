package it.cnr.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.workflow.Configuration;
import it.cnr.workflow.CryDetector;

public class TestCryDetector {

	public static void main(String[] args) throws Exception{
		
		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1.wav");
		
		
		CryDetector cryd = new CryDetector();
		List<File> outputAnnotations = cryd.detect(audio);
		cryd.eval(audio);
		
	}
	
}
