package it.cnr.evaluation.coroetal2024;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.sound.sampled.AudioFormat;

import it.cnr.clustering.DetectionManager;
import it.cnr.clustering.EnergyPitchFilterManager;
import it.cnr.features.CorpusCleaner;
import it.cnr.features.EnergyPitchFeatureExtractor;
import it.cnr.features.IslandDetector;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.workflow.configuration.AnomalousCryConfiguration;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class Evaluation {

	
	public static void evaluate(File[] inputFiles, File resultFile) throws Exception{
		
		
		
		AnnotatedReferenceFile result = new AnnotatedReferenceFile(resultFile, new File(resultFile.getAbsolutePath().replace(".wav", ".lab")));
		List<short[]> anomalousCry = result.extractSignals();
		int TP = 0;
		int FP = 0;
		int idd = 0;
		int nSamplesTP = 0;
		int nSamplesFP = 0;
		int maxI = 2;
		int nSamplesTP1 = 0;
		int nSamplesTP2 = 0;
		
		int ancrIdx = 0;
		for (short s[]:anomalousCry) {
			System.out.println("Analysing output annotation within "+Arrays.toString(result.annotationtime.get(idd)));
			boolean found = false;
			//int idx = 0;
			
			if (result.annotationtime.get(ancrIdx) [1]<10)
			{
				nSamplesTP1+=s.length;
				nSamplesTP+=s.length;
				TP++;
			}else if (result.annotationtime.get(ancrIdx) [1]<85) {
				nSamplesTP2+=s.length;
				nSamplesTP+=s.length;
				TP++;
			}else {
				nSamplesFP+=s.length;
				FP++;
			}
			
			/*
			for (File in:inputFiles) {
				System.out.println("SEARCHING FOR SOLUTION IN "+in.getName());
				File lab = new File(in.getAbsolutePath().replace(".wav", "_manual.lab"));
				AnnotatedReferenceFile arf = new AnnotatedReferenceFile(in, lab);
				String annotation = null;
				if (ancrIdx<6)
					annotation = "A";
				//String annotation = arf.search(s);
				if (annotation!=null) {
					System.out.println("Found annotation in "+in.getName()+" "+annotation);
					if (annotation.equals("A")) {
						TP++;
						nSamplesTP+=s.length;
					}else {
						FP++;
						nSamplesFP+=s.length;
					}
					found = true;	
					break;
				}else {
					System.out.println("No annotation associated - continuing...");
				}
				
				
				idx++;
				if(idx==1)
					break;
			}
			
			if (!found) {
				FP++;
				nSamplesFP+=s.length;
			}
			*/
			
			idd++;
			ancrIdx++;
		}
		
		System.out.println("FALSE NEGATIVE ASSESSMENT ... ");
		int id = 0;
		int FN = 0;
		int nSamplesFN = 0;
		for (File in:inputFiles) {
			System.out.println("Analysing file "+in.getName());
			
			File lab = new File(in.getAbsolutePath().replace(".wav", "_manual.lab"));
			AnnotatedReferenceFile arf = new AnnotatedReferenceFile(in, lab);
			List<short[]> anomalousCryF = arf.extractSignals();
			int TT = 0;
			for (int i=0;i<anomalousCryF.size();i++) {
				String ann = arf.annotations.get(i);
				
				if (ann.equals("A")) {
					TT+=anomalousCryF.get(i).length;
				}
			}
			
			if (id==0)
				nSamplesFN += TT-nSamplesTP1;
			else
				nSamplesFN += TT-nSamplesTP2;
			
			/*
			for (int i=0;i<anomalousCryF.size();i++) {
				String ann = arf.annotations.get(i);
				
				if (ann.equals("A")) {
					System.out.println("Searching for interval "+Arrays.toString(arf.annotationtime.get(i)));
					String found = result.search(anomalousCryF.get(i));
					if (found!=null && found.equals("A")) {
						System.out.println("Annotation ->"+found+" ");
					}else {
						System.out.println("False Negative found!");
						FN++;
						nSamplesFN+=anomalousCryF.get(i).length;
					}
				}
				
			}
			*/
			
			if(id==1)
				break;
			
			id++;
		}
		
		
		int summarySignalLength = result.signal.length;
		
		int totalSamples = 0;
		
		for (File in:inputFiles) {
			File lab = new File(in.getAbsolutePath().replace(".wav", "_manual.lab"));
			AnnotatedReferenceFile arf = new AnnotatedReferenceFile(in, lab);
			totalSamples+= arf.signal.length;
		}
		
		double segmentAccuracy = (double) (nSamplesTP)/ (double) (nSamplesTP+nSamplesFP);
		
		int nSamplesFPCD = summarySignalLength-nSamplesTP;
				
		double crydetectionAccuracy = (double) (totalSamples-nSamplesFPCD-nSamplesFN)/(double)totalSamples;
		
		double anomalousCrydetectionAccurcy = (double) (totalSamples-nSamplesFP-nSamplesFN)/ (double) totalSamples;//(double) (nSamplesTP+nSamplesTN)/(double)totalSamples;
		
		System.out.println("Segment accuracy (proportion of correctly identified samples): "+segmentAccuracy);
		System.out.println("Cry detection accuracy (proportion of anomalous samples captured by detection): "+crydetectionAccuracy);
		System.out.println("Actual accuracy (proportion of anomalous samples correctly selected): "+anomalousCrydetectionAccurcy);
		double precision = (double)nSamplesTP/(double) (nSamplesTP+nSamplesFP);
		double recall = (double)nSamplesTP/(double) (nSamplesTP+nSamplesFN);
		double f1 = 2*precision*recall/(precision+recall);
		
		System.out.println("Precision: "+precision);
		System.out.println("Recall: "+recall);
		System.out.println("F-measure: "+f1);
		
	}
	
	public static void main(String[] args) throws Exception{
		
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
				
				File ref = new File("./test_wave_files/islands_merged.wav");
				evaluate(allAudioToAnalyse, ref);
	}
		
	
}
