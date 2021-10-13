package it.cnr.tests;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import it.cnr.evaluation.Evaluator;

public class TestEvaluation {

	
	public static void main(String [] args) throws Exception {
		File audio = new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\TIN1.wav");
		
		List<File> waveFiles = new ArrayList<>();
		waveFiles.add(new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav"));
		waveFiles.add(new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_3.wav"));
		
		List<File> annotationFiles = new ArrayList<>();
		annotationFiles.add(new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_0.wav_energysegmentation_e0.1_p0.1\\audio_segment_0.wav_energysegmentation_e0.1_p0.1_LSTM_CTC.lab"));
		annotationFiles.add(new File("C:\\Users\\Utente\\Ricerca\\Experiments\\NINA Infant Cry\\TIN\\energysegmentation_0.5\\audio_segment_3.wav_energysegmentation_e0.1_p0.1\\audio_segment_3.wav_energysegmentation_e0.1_p0.1_LSTM_CTC.lab"));
		
		List<File> goldenAnnotations = new ArrayList<>();
		
		for (File w:waveFiles) {
			goldenAnnotations.add(new File(audio.getParentFile(),w.getName().replace(".wav", ".lab") ));
			
		}
		
		Evaluator evaluator = new Evaluator();
		evaluator.evaluate(annotationFiles, goldenAnnotations, waveFiles);
		
	}
	
	
}
