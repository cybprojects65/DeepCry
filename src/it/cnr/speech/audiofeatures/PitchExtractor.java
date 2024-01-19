package it.cnr.speech.audiofeatures;


import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import edu.cuny.qc.speech.AuToBI.core.Contour;
import edu.cuny.qc.speech.AuToBI.core.ContourIterator;
import edu.cuny.qc.speech.AuToBI.core.Pair;
import edu.cuny.qc.speech.AuToBI.core.WavData;
import edu.cuny.qc.speech.AuToBI.io.WavReader;
import it.cnr.speech.audiofeatures.AudioBits;

public class PitchExtractor {

	public double windowlength;
	public double windowshift;
	public int windowlengthsamples;
	public int windowshiftsamples;
	public int windowsamplesall;
	public int beginsamples;
	public int endsamples;
	public int spansamples;
	public double samples[];
	public Double pitchCurve[];
	
	public int samplerate;
	public String filepath;
	public boolean normalized;
	public double maxfrequency;
	public double maxfrequencywin;
	public double minfrequency;
	public boolean paintlines;
	public int waveheight;
	public double pitchWindowSec = 0.01d;
		
	public void setPitchWindowSec(double pitchWindowSec) {
		this.pitchWindowSec=pitchWindowSec; 
	}
	
	public static void main(String[] args) throws Exception{
		PitchExtractor extractor = new PitchExtractor();
		extractor.setPitchWindowSec(0.2);
		//System.out.println("is a question:"+extractor.isQuestionRough("./segments_v2\\cluster0\\audio_segment_37_word5.wav"));
		System.out.println("is a question:"+extractor.isQuestionRough("./segments_v2\\cluster0\\audio_segment_0_word0.wav"));
	}

	public boolean isQuestion(String wavefile) {
		calculatePitch(wavefile);
		
		if ((pitchCurve != null) && (pitchCurve.length > 0)) {
			// go backward until you find a non zero and non INF pitch p1
			int nonzeroidx = -1;
			for (int i = (pitchCurve.length - 1); i >= 0; i--) {
				if (pitchCurve[i]!=null && pitchCurve[i] != 0 && (pitchCurve[i]<1200) && !Double.isInfinite(pitchCurve[i]) && !Double.isNaN(pitchCurve[i])) {
					nonzeroidx = i;
					break;
				}
			}
			
			if (nonzeroidx == -1)
				return false;
			double pitchatidx=pitchCurve[nonzeroidx];
//			System.out.println("is question: found element at index "+nonzeroidx+" with value "+pitchatidx);
			// go backward until you find a pitch different from p1:p0
			int comparisonidx = -1;
			int tolerance = 5;
			for (int i = nonzeroidx-1; i >= 0; i--) {
				if ((pitchCurve[i]!=null && (pitchCurve[i] < pitchatidx-tolerance) || (pitchCurve[i] > pitchatidx+tolerance)) && 
						(pitchCurve[i]<1200) && !Double.isInfinite(pitchCurve[i]) && !Double.isNaN(pitchCurve[i])) {
					comparisonidx = i;
					break;
				}
			}	

			if (comparisonidx==-1)
				return false;
			double comparisonpitch=pitchCurve[comparisonidx];
//			System.out.println("is question: found comparison element at index "+comparisonidx+" with value "+comparisonpitch);
//			System.out.println(comparisonpitch+"->"+pitchatidx);
//			System.out.println("difference "+(comparisonpitch-pitchatidx));
			float k = 1.02f;
			//System.out.println("comparison "+comparisonpitch+" vs "+ pitchatidx +" "+(k*pitchatidx));
			//System.out.println("ratio "+comparisonpitch/pitchatidx+" vs 1.2");
			// if p1<p0 then it is a question
//			if (pitchatidx-comparisonpitch>0.7*pitchatidx)
//			if (comparisonpitch-pitchatidx>0.2*comparisonpitch)
//			if (comparisonpitch<pitchatidx){
			
			if ((comparisonpitch<k*pitchatidx)&&(comparisonpitch/pitchatidx<2)){
//				System.out.println("true");
				return true;
			}
			else{
//				System.out.println("false");
				return false;
			}
		}
		else return false;
		
	}

	public boolean isQuestionRough(String wavefile) {
		calculatePitch(wavefile);
		
		if ((pitchCurve != null) && (pitchCurve.length > 0)) {
			
			// go backward until you find a non zero and non INF pitch p1
						int nonzeroidx = -1;
						for (int i = (pitchCurve.length - 1); i >= 0; i--) {
							if (pitchCurve[i]!=null && pitchCurve[i] != 0 && (pitchCurve[i]<1200) && !Double.isInfinite(pitchCurve[i]) && !Double.isNaN(pitchCurve[i])) {
								nonzeroidx = i;
								break;
							}
						}
						
						if (nonzeroidx < 1)
							return false;
						
			Double pitchatidx=pitchCurve[nonzeroidx];
			Double comparisonpitch=pitchCurve[nonzeroidx-1];
			if (pitchatidx ==null || comparisonpitch == null) {
				if (pitchatidx!=null) 
					return true;
				else
					return false;
			}
			float k = 1.00f;
			//System.out.println("comparison "+comparisonpitch+" vs "+ pitchatidx +" "+(k*pitchatidx));
			//System.out.println("ratio "+comparisonpitch/pitchatidx+" vs 1.2");
			if ((comparisonpitch<k*pitchatidx)&&(comparisonpitch/pitchatidx<0.9)){
				return true;
			}
			else{
				return false;
			}
		}
		else return false;
		
	}
	
	public void getSamples(String filepath) {
		AudioBits audiosamplesreader = new AudioBits(new File(filepath));
		samples = audiosamplesreader.getDoubleVectorAudio();
		samplerate = (int)audiosamplesreader.getAudioFormat().getSampleRate();
		normalize(samples);
		audiosamplesreader.deallocateAudio();
		windowlength = 10;
		windowshift = 100;
	}

	
    public void calculatePitch(String wavefile) {
    	try{
    	WavReader reader = new WavReader();
    	WavData wave = reader.read(wavefile);
    	edu.cuny.qc.speech.AuToBI.PitchExtractor pe = new edu.cuny.qc.speech.AuToBI.PitchExtractor(wave);
    	Contour contour = pe.soundToPitch(pitchWindowSec,60,250);
    	int n = contour.size();
    	this.pitchCurve = new Double[n];
    	for (int i=0;i<n;i++) {
    		pitchCurve[i] = contour.get(i);
    		if (pitchCurve[i] == null || Double.isNaN(pitchCurve[i]) || Double.isInfinite(pitchCurve[i]))
    			pitchCurve[i] = 0d;
    	}
    	
    	
    	/*
    	ContourIterator iterator = contour.iterator();
    	int i = 0;
    	while (iterator.hasNext()){
    		Pair<Double,Double> pair = iterator.next();
    		
    		if (pair.second == null){
    			System.out.println(pair.second);
    			pair.second=Double.NaN;
    		}
    		pitchCurve[i]=pair.second;
    		i++;
    	}
    	System.out.println("I:"+i);
    	*/
    	
    	}catch(Exception e ){}
    }
    

    public void calculatePitch(String wavefile, double minPitch, double maxPitch) {
    	try{
    	WavReader reader = new WavReader();
    	WavData wave = reader.read(wavefile);
    	double duration = wave.getDuration();
    	
    	edu.cuny.qc.speech.AuToBI.PitchExtractor pe = new edu.cuny.qc.speech.AuToBI.PitchExtractor(wave);
    	Contour contour = pe.soundToPitch(pitchWindowSec,minPitch,maxPitch);
    	 //pe.soundToPitchAc(pitchWindowSec, minPitch, 1.0, 15, 0.03, 0.45, 0.01, 0.35, 0.14, maxPitch);
    	int contentSize = contour.size();
    	//HashMap<Integer,Double> pitchMap = new HashMap<>();
    	int pitchCurveLen = (int) (duration/pitchWindowSec);
    	
    	this.pitchCurve = new Double[pitchCurveLen];
    	
    	for (int i = 0; i < contentSize; ++i) {
    		double time = contour.timeFromIndex(i)-pitchWindowSec; // the indicated time is the final marker of the pitch
    		double value = contour.get(i);
            //System.out.println("point[" + i + "]: " + contour.get(i) + " -- " + contour.timeFromIndex(i));
            
            Integer index = (int) (time/pitchWindowSec);
            pitchCurve[index] = value;
          }
    	
    	System.out
        .println(
            "pitch elements "+pitchCurve.length);
    	
    	/*
    	double time = 0;
    	List<Double> pitchArray = new ArrayList<>();
    	int i = 0;
    	while (i<contentSize)
    	{
    		double value = contour.get(time);
    		contour.
    		System.out.println("time "+time+ " value " + value);
    		
    		
    		pitchArray.add(value);
    		time = time + pitchWindowSec;
    		i++;
    	}
    	
    	this.pitchCurve = new Double[pitchArray.size()];
    	pitchCurve = pitchArray.toArray(pitchCurve);
    	*/
    	/*
    	this.pitchCurve = new Double[contour.size()];
    	ContourIterator iterator = contour.iterator();
    	int i = 0;
    	while (iterator.hasNext()){
    		Pair<Double,Double> pair = iterator.next();
    		System.out.println("ptc "+pair.first+" "+pair.second);
    		if (pair.second == null){
    			System.out.println(pair.second);
    			pair.second=Double.NaN;
    		}
    		pitchCurve[i]=pair.second;
    		i++;
    	}
    	
    	*/
    	
    	
    	}catch(Exception e ){
    		e.printStackTrace();
    		
    	}
    }
    
	public double[] normalize(double af[]) {
		double f = 0.0F;
		for (int i = 0; i < af.length; i++) {
			if (f < Math.abs(af[i])) {
				f = Math.abs(af[i]);
			}
		}

		for (int j = 0; j < af.length; j++) {
			af[j] *= 128D / (double) f;
		}

		normalized = true;
		return af;
	}

	public static float[] correlate(int size, float[] x){
		    float[] R = new float[size];
		    float sum;
		    float af2[] = hammingWindow(size);
		    for (int i=0;i<size;i++) {
		        sum=0;
		        for (int j=0;j<size-i;j++) {
		            sum+=x[j]*x[j+i]*af2[j];
		        }
		        R[i]=sum;
		    }
		    
		    return R;
		}

	public float[] autoCorrelate(float af[], int i, int j) {
		int k = af.length;
		if (k < i + j) {
			return null;
		}
		float af1[] = new float[i];
		float af2[] = hammingWindow(j);
		for (int l = 0; l < i; l++) {
			for (int i1 = 0; i1 < j; i1++) {
				af1[l] += af[i1] * af[i1 + l] * af2[l];
			}

		}

		return af1;
	}

	public static float[] hammingWindow(int i) {
		float af[] = new float[i];
		for (int j = 0; j < i; j++) {
			af[j] = 0.5F - 0.5F * (float) Math.cos((6.2831853071795862D * (double) j) / (double) i);
		}

		return af;
	}

	public void outPitches() {
		System.out.println(pitchCurve.length + " Pitches");
		for (int i = 0; i < pitchCurve.length; i++) {
			System.out.print(pitchCurve[i] + ", ");
		}

	}

	public double[] smoothWithDegreeFour(double af[]) {
		int i = af.length;
		af[3] = (af[0] + af[1] + af[2] + af[3] + af[4]) / 5F;
		af[2] = (af[0] + af[1] + af[2] + af[3]) / 4F;
		af[1] = (af[0] + af[1] + af[2]) / 3F;
		af[0] = (af[0] + af[1]) / 2.0F;
		af[i - 4] = (af[i - 1] + af[i - 2] + af[i - 3] + af[i - 4] + af[i - 5]) / 5F;
		af[i - 3] = (af[i - 1] + af[i - 2] + af[i - 3] + af[i - 4]) / 4F;
		af[i - 2] = (af[i - 1] + af[i - 2] + af[i - 3]) / 3F;
		af[i - 1] = (af[i - 1] + af[i - 2]) / 2.0F;
		for (int j = 4; j < i - 4; j++) {
			af[j] = (af[j] + af[j + 1] + af[j - 1] + af[j + 2] + af[j - 2] + af[j + 3] + af[j - 3] + af[j + 4] + af[j - 4]) / 9F;
		}

		return af;
	}

}
