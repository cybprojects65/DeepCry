package it.cnr.speech.filters;

public class HighPassFilter {
	
	double samplingRate; double cutoffFrequency;
	
    // Constructor
    public HighPassFilter(double samplingRate, double cutoffFrequency) {
    	this.samplingRate = samplingRate;
    	this.cutoffFrequency = cutoffFrequency;
    }
    
    public double[] highPass(double [] inputSignal) {
    	LowPassFilterDynamic lowpass = new LowPassFilterDynamic(cutoffFrequency, samplingRate, 0.01);

    	double[] lowPassFilteredSignal =  lowpass.applyFilter(inputSignal);
    	
    	// Subtract the low-pass filtered signal from the original signal to get the high-pass filtered signal
        int signalLength = inputSignal.length;
        double[] highPassFilteredSignal = new double[signalLength];
        for (int i = 0; i < signalLength; i++) {
            highPassFilteredSignal[i] = inputSignal[i] - lowPassFilteredSignal[i];
        }

        return highPassFilteredSignal;

    }
    
    
}