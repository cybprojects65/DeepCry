package it.cnr.tests.coroetal2023;

import java.io.File;

import it.cnr.workflow.coroetal2023.CryDetectorOptimiserCoroetAl2023;

public class TestCryDetectorSearcher {

	public static void main(String[] args) throws Exception{
		//Logger rootLogger = (Logger)LoggerFactory.getLogger(Logger.ROOT_LOGGER_NAME);
		
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1\\TIN1.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN2\\TIN2.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB1\\SUB1.wav");
		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB2\\SUB2.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB3\\SUB3.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB4\\SUB4.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\SUB\\SUB5\\SUB5.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\NIDO\\NIDO1\\NIDO1.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\NIDO\\NIDO2\\NIDO2.wav");
		//File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\NIDO\\NIDO3\\NIDO3.wav");
		CryDetectorOptimiserCoroetAl2023 cryd = new CryDetectorOptimiserCoroetAl2023();
		cryd.findOptimalConfiguration(audio);
		
	
	}
	
}
