package it.cnr.speech.filters;

import org.apache.commons.math3.complex.Complex;

public class GreenwoodFilterBank extends LowPassFilterDynamic{
	
	public final static double lowerFilterFreq = 165.4;//200; //minimum Greenwood filterbank frequency
	public int numMelFilters = 15;//17; //at 13 this will be under 4000 to simulate Greenwood frequencies
	int melFreqIndex;
	double highestFrequency;
	
	public GreenwoodFilterBank(int melFreqIndex, double samplingRate, double windowSizeSec, int numFilters, double highestFrequency) {
		super(0, samplingRate, windowSizeSec);
		this.numMelFilters = numFilters;
		this.melFreqIndex = melFreqIndex;
		this.highestFrequency = highestFrequency;
	}
	    
    protected static double log10(double value){
        return Math.log(value) / Math.log(10);
    }
    private static double inverseMel(double x){
        double temp = Math.pow(10, x / 2595) - 1;
        return 700 * (temp);
    }
    protected static double freqToMel(double freq){
        return 2595 * log10(1 + freq / 700.0);
        
    }

    protected static double freqToGF(double freq){
        return 165.4 * (1 + (freq / 634.5));
        
    }
    
    private static double inverseGF(double gf){
        return 634.5* ( (gf/165.4)-1);
    }
    
    public static double greenwoodFunction(double frequency, double a, double fMin) {
        return a * Math.log10(1 + frequency / fMin);
    }
    
    private double centerFreq(int i,double samplingRate){
    	
    	
        double mel[] = new double[2];
        mel[0] = freqToMel(lowerFilterFreq);
        //mel[1] = freqToMel(samplingRate / 2);
        mel[1] = freqToMel(highestFrequency);
        // take inverse mel of:
        double temp = mel[0] + ((mel[1] - mel[0]) / (numMelFilters + 1)) * i;
        return inverseMel(temp);
    	
    	
    	/*
    	double gf0 = freqToGF(lowerFilterFreq);
    	double gf1 = freqToGF(highestFrequency);
    	double step = (gf1-gf0)/ (numMelFilters+1);
    	
    	double gfreq = gf0+( (i*step)+step);
    	return inverseGF(gfreq);
    	*/
    	/*
    	double startFrequency = 20.0; // Start from 20 Hz
        double endFrequency = 8000.0; // End at 8000 Hz

        // Set parameters for the Greenwood function
        double a = 165.4; // Constant for Greenwood function
        double fMin = 165.4; // Minimum frequency mapped to the cochlea
        
        double frequency = startFrequency + (endFrequency - startFrequency) * i / (numMelFilters - 1);
        double cochleaPositions = greenwoodFunction(frequency, a, fMin);
        return cochleaPositions;
        */
    }
    
    
    private int calcCbin(int melFreqIndex) {
    
    	double fc = centerFreq(melFreqIndex,samplingRate);
    	
        int cbin = (int) Math.floor(fc * (double) windowSize / samplingRate);
        
        
    	return cbin;
    }
    
    private boolean displayed = false; 
    protected Complex[] applyBandpassFilter(Complex[] complexSpectrum) {
    	
    	
    	int cbin = calcCbin(melFreqIndex);
    	int cbinprev = calcCbin(melFreqIndex-1);
    	int cbinpost = calcCbin(melFreqIndex+1);
    	
    	if (!displayed) {
    		System.out.println("Filterbank index: "+melFreqIndex);
    		System.out.println("Corresp freq: "+centerFreq(melFreqIndex,samplingRate));
    		System.out.println("Corresp index in FFT: "+cbin);
    		
    		System.out.println("Previous freq: "+centerFreq(melFreqIndex-1,samplingRate));
    		System.out.println("Corresp index in FFT: "+cbinprev);
    		
    		System.out.println("Next freq: "+centerFreq(melFreqIndex+1,samplingRate));
    		System.out.println("Corresp index in FFT: "+cbinpost);
    		
    		System.out.println("Filterbank indexes: ["+cbinprev+" ; "+cbin+" ; "+cbinpost+"]");
    		
    		displayed = true;
    	}
    	
    	
        Complex[] filteredSpectrum = new Complex[complexSpectrum.length];
        
        for (int i = 0; i < cbinprev; i++)
        	filteredSpectrum[i] = Complex.ZERO;
        		
        for (int i = cbinprev; i <= cbin; i++){
            double weight = ((i - cbinprev + 1) / (cbin - cbinprev + 1));
        	filteredSpectrum[i] = complexSpectrum[i].multiply(weight);
        }
        
        for (int i = cbin + 1; i <= cbinpost; i++){
        	double weight = (1 - ((i - cbin) / (cbinpost - cbin + 1)));
        	filteredSpectrum[i] = complexSpectrum[i].multiply(weight);
        }
        
        for (int i = (cbinpost+1); i < filteredSpectrum.length; i++)
        	filteredSpectrum[i] = Complex.ZERO;
        
        return filteredSpectrum;
    }

    
    public static void main(String[] args) {
    	
    	int numMelFilters = 8;
    	double samplingFrequency = 16000;
    	double maxFrequency = 3000;//8000;
    	
    	for (int i= 0; i<numMelFilters;i++){
    		
    		GreenwoodFilterBank mfb = new GreenwoodFilterBank(i, 16000.0, 0.01, numMelFilters, maxFrequency);
    		int cbin = mfb.calcCbin(i);
    		double fcu = mfb.centerFreq(i, 16000);
    		
    		int cbinprev = mfb.calcCbin(i-1);
    		double fpr = mfb.centerFreq(i-1, 16000);
    		int cbinpost = mfb.calcCbin(i+1);
    		double fpo = mfb.centerFreq(i+1, 16000);
    		
    		System.out.println("MF N.: "+i+" ("+fpr+" ; "+fcu+" ; "+fpo+")"+" -> ["+cbinprev+" ; "+cbin+" ; "+cbinpost+"]");
    		System.out.println("Freq: "+mfb.centerFreq(i,16000));
    	}
    }
    
    
}

