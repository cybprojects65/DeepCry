package it.cnr.speech.filters;

import java.util.Arrays;

import org.apache.commons.math3.complex.Complex;
import org.apache.commons.math3.transform.DftNormalization;
import org.apache.commons.math3.transform.FastFourierTransformer;

import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;

public class LowPassFilterDynamic {
	
	protected final static double lowerFilterFreq = 200; //minimum Greenwood frequency = lowest audible frequency
	
	double freq; double samplingRate; int windowSize;  int cbin; 
	
	
	public LowPassFilterDynamic(double freq, double samplingRate, double windowSizeSec) {
		this.freq = freq;
		this.samplingRate = samplingRate;
		int windowSizeSamples = SignalProcessing.timeToSamples(windowSizeSec, samplingRate);
		this.windowSize = UtilsMath.powerTwoApproximation(windowSizeSamples);
		this.cbin = (int) Math.floor(freq * (double) windowSize / samplingRate);
		System.out.println("Window size to be used in filter: "+this.windowSize+" (samples)");
	}
	
    public double[] applyFilter(double[] signal) {
    	
    	//double windowSizeSec = 0.01; //ms
    	int signalLength = signal.length;
        int hopSize = windowSize / 2;

        double[] filteredSignal = new double[signalLength];

        for (int i = 0; i < signalLength; i += hopSize) {
        	if ((i+windowSize)>signalLength)
        		break;
            // Extract a windowed segment of the signal
            double[] windowedSegment = getWindowedSegment(signal, i, windowSize);
            windowedSegment = hammingWindow(windowedSegment);
            
            // Compute the Fourier Transform for the windowed segment
            Complex[] complexSpectrum = computeFourierTransform(windowedSegment);

            // Apply bandpass filter to the spectrum
            Complex[] filteredSpectrum = applyBandpassFilter(complexSpectrum);
            
            // Reconstruct the windowed segment using Inverse Fourier Transform
            double[] reconstructedSegment = reconstructSignal(filteredSpectrum);

            // Overlap and add the reconstructed segment to the filtered signal
            overlapAndAdd(filteredSignal, reconstructedSegment, i);
            
            
        }

        return filteredSignal;
    }

    public static double[] hammingWindow(double[] windowedSegment){
        double w[] = new double[windowedSegment.length];
        for (int n = 0; n < windowedSegment.length; n++){
            w[n] = 0.54 - 0.46 * Math.cos( (2 * Math.PI * n) / (double) (windowedSegment.length - 1) );
        }

        double[] windowedSegmentHW = new double[windowedSegment.length];
        
            for (int n = 0; n < windowedSegment.length; n++){
            	windowedSegmentHW[n] = windowedSegment[n]* w[n];
            }
        return windowedSegmentHW;
    }
    
    public static double[] getWindowedSegment(double[] signal, int startIdx, int windowSize) {
        int signalLength = signal.length;
        //int endIdx = Math.min(startIdx + windowSize, (signalLength-1));

        double[] windowedSegment = new double [windowSize];
        
        for (int i =0 ;i<windowSize;i++) {
        	int nextIdx = startIdx+i;
        	if (nextIdx< (signal.length-1 ) )
        		windowedSegment [i] = signal [nextIdx];
        	else
        		windowedSegment [i] = 0;
        }
        //double[] windowedSegment = Arrays.copyOfRange(signal, startIdx, endIdx);

        // Apply window function (e.g., Hamming window) if needed

        return windowedSegment;
    }

    public static Complex[] computeFourierTransform(double[] signal) {
        FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        return transformer.transform(signal, org.apache.commons.math3.transform.TransformType.FORWARD);
    }

    protected Complex[] applyBandpassFilter(Complex[] complexSpectrum) {
    	
        Complex[] filteredSpectrum = new Complex[complexSpectrum.length];
        
        for (int i = 0; i <= cbin; i++)
        	filteredSpectrum[i] = complexSpectrum[i];
        
        for (int i = (cbin+1); i < complexSpectrum.length; i++)
        	filteredSpectrum[i] = Complex.ZERO;
        		
        return filteredSpectrum;
    }

    protected double[] reconstructSignal(Complex[] complexSpectrum) {
        FastFourierTransformer transformer = new FastFourierTransformer(DftNormalization.STANDARD);
        Complex[] inverseTransform = transformer.transform(complexSpectrum, org.apache.commons.math3.transform.TransformType.INVERSE);

        // Extract the real part of the inverse transform as the reconstructed signal
        return Arrays.stream(inverseTransform).mapToDouble(Complex::getReal).toArray();
    }

    protected void overlapAndAdd(double[] outputSignal, double[] segment, int startIdx) {
        int segmentLength = segment.length;
        int outputLength = outputSignal.length;

        for (int i = 0; i < segmentLength; i++) {
            if (startIdx + i < outputLength) {
                outputSignal[startIdx + i] += segment[i];
            }
        }
    }

   
}

