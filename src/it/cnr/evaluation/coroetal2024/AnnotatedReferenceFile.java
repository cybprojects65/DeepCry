package it.cnr.evaluation.coroetal2024;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.workflow.utils.SignalProcessing;
import jdk.internal.misc.Signal;

public class AnnotatedReferenceFile {

	double samplingFrequency;
	short [] signal;
	List<double []> annotationtime;
	List<String> annotations;
	
	public AnnotatedReferenceFile(File audioFile, File annotation) throws Exception{
		
		AudioBits ab = new AudioBits(audioFile);
		samplingFrequency = ab.getAudioFormat().getSampleRate();
		signal = ab.getShortVectorAudio();
		ab.ais.close();
		
		parseAnnotation(annotation);
	}
	
	public List<short[]> extractSignals(){
		
		List<short[]> signalsExtracted = new ArrayList<>(); 
		
		for (double[] annotationT:annotationtime) {
			double t0 = annotationT[0];
			double t1 = annotationT[1];
			int s0 = SignalProcessing.timeToSamples(t0, samplingFrequency);
			int s1 = SignalProcessing.timeToSamples(t1, samplingFrequency);
			short[] sig = Arrays.copyOfRange(signal, s0,s1+1);
			signalsExtracted.add(sig);
		}
		
		return signalsExtracted;
		
	}
	
	
	public void parseAnnotation(File annotationFile) throws Exception{
		List<String> annotationsF = Files.readAllLines(annotationFile.toPath());
		annotationtime = new ArrayList<double[]>();
		annotations = new ArrayList<String>();
		
		for (String annotation:annotationsF) {
			
			String elements[] = annotation.split(" ");
			double t0 = Double.parseDouble(elements[0]);
			double t1 = Double.parseDouble(elements[1]);
			double [] t = {t0, t1};
			String an = "";
			if (elements.length>2)
				an = elements[2];
			if (an.trim().length()>0) {
				annotationtime.add(t);
				annotations.add(an);
			}
		}
	}
	
	//if the interval intersect a known interval for more than 30% then return the annotation 
	public String getAnnotation(double t0,double t1){
		
		int idx = 0;
		for (double [] interval:annotationtime) {
			double i0 = interval[0];
			double i1 = interval[1];
			double intersection = 0;
			if ( t0>i0 && t0<i1 && t1>i1)
				intersection = i1-t0;
			else if ( t0>i0 && t0<i1 && t1<i1)
				intersection = t1-t0;
			else if ( t0<i0 && t1>i0 && t1<i1)
				intersection = t1-i0;
			else if ( t0<=i0 && t1>=i1)
				intersection = i1-i0;
			
			if (intersection>0) {
				double perc = (intersection/(t1-t0));
				if (perc>0.2) {
					return annotations.get(idx);
				}
			}
			idx++;
		}
		
		
		return null;
	}
	
	public static short[] removezeros(short[] signal) {
		
		List<Short> s = new ArrayList<>();
		
		for (short d: signal) {
			if (d!=0)
				s.add(d);
		}
		
		short red []=new short[s.size()];
		int i=0;
		for (Short d:s) {
			red[i] = d;
			i++;
		}
		
		return red;
	} 
	
	public static boolean equals(short []s1,short s2[]) {
		
		short [] sr1 = s1; //removezeros(s1);
		short [] sr2 = s2; //removezeros(s2);
		/*
		if (sr1.length!=sr2.length) {
			return false;
		}
		*/
		int equal = 0;
		for (int i=0;i<sr1.length;i++) {
			if (sr1[i] == sr2[i]) {
				equal ++;
			}
		}
		double score = (double) equal/ (double) sr1.length;
		if (score>0.01)
			return true;
		else
			return false;
	}
	
	public String search(short []s) {
		
		List<short[]> extracted = extractSignals();
		int id=0;
		for (short[] se:extracted) {
				
			short s1 []= s;
			short s2 []= se;
			if (s.length<se.length) {
				s2 = s;
				s1 = se;
			}
			
			for (int i=0;i<s1.length;i++) {
				if ((i+s2.length)<s1.length) {
					short ss [] = Arrays.copyOfRange(s1, i, (i+s2.length));
					if (equals(ss,s2)) {
						double t [] = annotationtime.get(id);
						System.out.println("Found evidence between ["+t[0]+";"+t[1]+"]");
						return getAnnotation(t[0], t[1]);
					}
				}
			}
			
			id++;
		}
		
		return null;
		
	}

	public String search1(short []s) {
		
		double t0=0;
		double t1=0;
		//s = removezeros(s);
		for (int i=0;i<signal.length;i++) {
			t0 = SignalProcessing.samplesToTime(i, samplingFrequency);
			t1 = SignalProcessing.samplesToTime((i+s.length), samplingFrequency);
			if ((i+s.length)<signal.length) {
				short ss [] = Arrays.copyOfRange(signal, i, (i+s.length));
				//ss = removezeros(ss);
				if (equals(ss,s)) {
					System.out.println("Found evidence between ["+t0+";"+t1+"]");
					return getAnnotation(t0, t1);
				}
			}
		}
		
		return null;
		
	}
	
	
	//if the interval intersect a known interval for more than 30% then return the annotation 
	public String searchAnnotation(double t0,double t1){
			
			int idx = 0;
			for (double [] interval:annotationtime) {
				double i0 = interval[0];
				double i1 = interval[0];
				double intersection = 0;
				if ( t0>i0 && t0<i1 && t1>i1)
					intersection = i1-t0;
				else if ( t0>i0 && t0<i1 && t1<i1)
					intersection = t1-t0;
				else if ( t0<i0 && t1>i0 && t1<i1)
					intersection = t1-i0;
				else if ( t0<=i0 && t1>=i1)
					intersection = i1-i0;
				
				if (intersection>0) {
					double perc = (intersection/(t1-t0));
					if (perc>0.2) {
						return annotations.get(idx);
					}
				}
				idx++;
			}
			
			
			return null;
		}
		
}
