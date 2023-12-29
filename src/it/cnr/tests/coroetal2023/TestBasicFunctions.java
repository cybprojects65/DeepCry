package it.cnr.tests.coroetal2023;

import java.io.File;

import it.cnr.evaluation.Evaluator;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.Utils;

public class TestBasicFunctions {

	
	public static void main(String[] args) throws Exception {
		
		double [] e1 = {1,-1};
		//double [] e2 = {0,1};
		double [] e2 = {1,0};
		double angle = Utils.angle(e1, e2);
		
		System.out.println("Angle="+angle);
		
		System.out.println("K ALL="+Evaluator.CohensKappa(7521595, 14050791, 1113881, 6259671));
		System.out.println("K NUR="+Evaluator.CohensKappa(2989868, 3522073, 2518222, 362095));
		System.out.println("K SUB="+Evaluator.CohensKappa(3972846, 5538998, 1671965, 412106));
		System.out.println("K INT="+Evaluator.CohensKappa(757659, 4035718, 1337095, 129948));
		
		
		File o1 = new File ("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\clustering.lab");
		File o2 = new File ("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\clustering_uni.lab");
		CorpusCleaner.uniformLabelSegments(o1, o2);
		
		o1 = new File ("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_3.wav_energysegmentation_e0.1_p0.1\\LSTM.lab");
		o2 = new File ("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_3.wav_energysegmentation_e0.1_p0.1\\LSTM_uni.lab");
		
		CorpusCleaner.uniformLabelSegments(o1, o2);
		
		
	}
	
}
