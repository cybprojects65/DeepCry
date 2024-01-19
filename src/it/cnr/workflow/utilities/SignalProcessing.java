package it.cnr.workflow.utilities;

import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import javax.swing.JPanel;

import org.apache.commons.math3.complex.Complex;

import com.rapidminer.example.Attribute;
import com.rapidminer.example.Example;
import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.MemoryExampleTable;

import it.cnr.clustering.BigSamplesTable;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.filters.Delta;
import it.cnr.speech.filters.LowPassFilterDynamic;
import marytts.signalproc.display.SpectrogramCustom;
import marytts.signalproc.window.Window;

/**
 * includes tools for basic signal transformations: delta + double delta center frequency cepstral coefficients calculation spectrum frequency cut transformation to and from Rapid Miner Example Set filterbanks fequency to mel frequency to index in fft sinusoid signal generation inverse mel log10 mel filterbanks sample to time and time to sample signal timeline generation index to time in spectrogram spectrogram calculation and display time to index in spectrogram
 * 
 * @author coro
 * 
 */
public class SignalProcessing {

	public static double calculateSpectralEntropy(double[] frequencies) {
        // Normalize the spectrum to obtain probabilities
        double[] probabilities = normalizeSpectrum(frequencies);

        // Calculate entropy
        double entropy = 0.0;
        int i=0;
        for (double probability : probabilities) {
            if (probability > 0.0) {
                entropy -= probability * Math.log(probability);
            	//double invdB = Math.pow(10, frequencies[i]/10d); 
            	//entropy += invdB;
            }
            i++;
        }

        return entropy;
    }
	
	 public static double[] normalizeSpectrum(double[] frequencies) {
		 
		 	//double min = Utils.getMin(frequencies);
	        // Calculate the sum of frequencies
	        double sum = 0.0;
	        for (double frequency : frequencies) {
	        	double invdB = Math.pow(10, frequency/10d); 
	            //sum += frequency;//frequency+min+1;
	        	sum += invdB;
	        }

	        // Normalize frequencies to obtain probabilities
	        double[] probabilities = new double[frequencies.length];
	        for (int i = 0; i < frequencies.length; i++) {
	        	double invdB = Math.pow(10, frequencies[i]/10d);
	            probabilities[i] = invdB / sum;//frequencies[i] / sum;//Math.abs(frequencies[i]) / sum;
	        }

	        return probabilities;
	    }
	public static List<double[]> getAnnotatedSignalSegments(File audio, double times[], String[] labels) throws Exception{
		
		AudioBits ab = new AudioBits(audio);
		double[] signal = ab.getDoubleVectorAudio();
		int sfrequency = (int) ab.getAudioFormat().getSampleRate();
		ab.ais.close();
		List<double[]> segments = new ArrayList();
		int totalSamples=0;
		
		for (int i=0;i<times.length;i++) {
			String label = labels[i].trim();
			if (label.length()>0) {
				
				double time0 = times[i];
				double time1 = times[times.length-1];
				if (i< (times.length-1 ) )
					time1 = times[i+1];
				
				int sampletime0 = timeToSamples(time0, (double) sfrequency);
				int sampletime1 = timeToSamples(time1, (double) sfrequency);
				
				double[] subsignal = Arrays.copyOfRange(signal, sampletime0, sampletime1+1);
				totalSamples=totalSamples+subsignal.length;
				
				segments.add(subsignal);
			}
		}
		
		
		
		return segments;
		
	}

	public static short[] extractAnnotatedSignal(File audio, double times[], String[] labels) throws Exception{
		
		AudioBits ab = new AudioBits(audio);
		short[] signal = ab.getShortVectorAudio();
		int sfrequency = (int) ab.getAudioFormat().getSampleRate();
		ab.ais.close();
		List<short[]> segments = new ArrayList();
		int totalSamples=0;
		for (int i=0;i<times.length;i++) {
			String label = labels[i].trim();
			if (label.length()>0) {
				
				double time0 = times[i];
				double time1 = times[times.length-1];
				if (i< (times.length-1 ) )
					time1 = times[i+1];
				
				int sampletime0 = timeToSamples(time0, (double) sfrequency);
				int sampletime1 = timeToSamples(time1, (double) sfrequency);
				
				short[] subsignal = Arrays.copyOfRange(signal, sampletime0, sampletime1+1);
				totalSamples=totalSamples+subsignal.length;
				
				segments.add(subsignal);
			}
		}
		
		
		double reduction_perc = 100d*(double)(signal.length-totalSamples)/(double)signal.length;
		System.out.println("The original signal will be reduced of "+reduction_perc+"%");
		short[] unifiedSignal = new short[totalSamples];

		int idx = 0;
		for (int i = 0;i<segments.size();i++) {
				
			short[] singlesignal = segments.get(i);
			for (int j=0;j<singlesignal.length;j++) {
				unifiedSignal[idx] = singlesignal[j];
				idx = idx+1;
			}
		}
		
		return unifiedSignal;
		
	}
	
	
	public static double calculateAverageEnvelopeLevel(double[] signal) {
        double sum = 0.0;

        for (double sample : signal) {
            sum += Math.abs(sample);
        }

        return sum / signal.length;
    }
	
	public static double[] normaliseByAverageEnvelopeLevel (double signal[]) {
		
		double ael = calculateAverageEnvelopeLevel(signal);
		double aelSignal [] = new double[signal.length];
		for (int i=0;i<signal.length;i++) {
			aelSignal[i] = signal[i]/ael;
		}
		return aelSignal;
	}
	
	public int windowShiftSamples;
	public int windowSizeSamples;
	public int samplingRate;
	public double signal[];
	public void getSignal(File audio) throws Exception{
		AudioBits bits = new AudioBits(audio);
		 signal = bits.getDoubleVectorAudio();
		samplingRate = (int) bits.getAudioFormat().getSampleRate(); 
		bits.ais.close();
	}
	
	public double[][] shortTermFFT(File audio, double windowSize, double windowShift) throws Exception{
		
		getSignal(audio);
		return shortTermFFT(signal, samplingRate ,windowSize, windowShift);
	}
	
	public double[][] shortTermFFT(double[] signal, int samplingRate ,double windowSize, double windowShift) throws Exception{	
		windowSizeSamples = SignalProcessing.timeToSamples(windowSize, samplingRate);
		System.out.println("Original window sample: "+windowSize+"s"+" "+windowSizeSamples+" (samples)");
		windowSizeSamples = UtilsMath.powerTwoApproximation(windowSizeSamples);
		System.out.println("Approx window sample: "+SignalProcessing.samplesToTime(windowSizeSamples,samplingRate) +"s"+" "+windowSizeSamples+" (samples)");
		
		windowShiftSamples = SignalProcessing.timeToSamples(windowShift, samplingRate);
		windowShiftSamples = UtilsMath.powerTwoApproximation(windowShiftSamples);
		
		System.out.println("Running FFT with "+windowSizeSamples+" by "+windowShiftSamples+" ...");
		
		List<double[]> spectra = new ArrayList<double[]>();
        
        for (int i = 0; i < signal.length; i += windowShiftSamples) {
        	if ((i+windowSizeSamples)>signal.length)
        		break;
        	
            // Extract a windowed segment of the signal
            double[] windowedSegment = LowPassFilterDynamic.getWindowedSegment(signal, i, windowSizeSamples);
            windowedSegment = LowPassFilterDynamic.hammingWindow(windowedSegment);
            
            // Compute the Fourier Transform for the windowed segment
            Complex[] complexSpectrum = LowPassFilterDynamic.computeFourierTransform(windowedSegment);

            // Apply bandpass filter to the spectrum
            double[] absSpectrum = new double[complexSpectrum.length];
            for (int k=0;k<absSpectrum.length;k++) {
            	
            	absSpectrum[k] = complexSpectrum[k].abs();
            }
            
            spectra.add(absSpectrum);
        }

        double[][] spectrum = spectra.toArray(new double[spectra.size()][]);
		/*
        ShortTermSpectrumAnalyser spectrumAnalyser = new ShortTermSpectrumAnalyser
                (new BufferedDoubleDataSource(signal), windowSizeSamples, Window.get(Window.HAMMING, windowSizeSamples), windowShiftSamples, samplingRate);
        System.out.println("Initialised");
        List<double[]> spectra = new ArrayList<double[]>();
        long startTime = System.currentTimeMillis();
        FrameBasedAnalyser.FrameAnalysisResult<double[]>[] results = spectrumAnalyser.analyseAllFrames();
        long endTime0 = System.currentTimeMillis();
        System.out.println("Computed " + spectra.size() + " spectra in " + (endTime0-startTime) + " ms.");
        
        for (int i=0; i<results.length; i++) {
            double[] spectrum = (double[]) results[i].get();
            System.out.println("sp length "+spectrum.length);
            spectra.add(spectrum);
        }
        long endTime1 = System.currentTimeMillis();
        System.out.println("Added " + spectra.size() + " spectra in " + (endTime1-startTime) + " ms.");
        
        double[][] spectrum = spectra.toArray(new double[spectra.size()][]);
        */
        
        
        
        return spectrum;
        /*
            // Frequency resolution of the FFT:
            deltaF = spectrumAnalyser.getFrequencyResolution(); 
            long startTime = System.currentTimeMillis();
            spectra_max = Double.NaN;
            spectra_min = Double.NaN;
            FrameBasedAnalyser.FrameAnalysisResult<double[]>[] results = spectrumAnalyser.analyseAllFrames();
            for (int i=0; i<results.length; i++) {
                double[] spectrum = (double[]) results[i].get();
                spectra.add(spectrum);
                // Still do the preemphasis inline:
                for (int j=0; j<spectrum.length; j++) {
                    double freqPreemphasis = PREEMPHASIS / Math.log(2) * Math.log((j+1)*deltaF/1000.);
                    spectrum[j] += freqPreemphasis;
                    if (Double.isNaN(spectra_min) || spectrum[j] < spectra_min) {
                        spectra_min = spectrum[j];
                    }
                    if (Double.isNaN(spectra_max) || spectrum[j] > spectra_max) {
                        spectra_max = spectrum[j];
                    }
                }
            }
            long endTime = System.currentTimeMillis();
            System.err.println("Computed " + spectra.size() + " spectra in " + (endTime-startTime) + " ms.");

            spectra_indexmax = (int) (FREQ_MAX / deltaF);
            if (spectra_indexmax > fftSize/2)
                spectra_indexmax = fftSize/2; // == spectra[i].length
                */
	}
	
	public static double[][] addDeltaDouble(double[][] features) throws Exception {
		int vectorL = features[0].length;
		double[][] delta = new double[features.length][features[0].length * 3];

		for (int k = 0; k < features.length; k++) {
			for (int g = 0; g < vectorL; g++) {
				delta[k][g] = features[k][g];
			}
		}

		Delta.calcDelta(delta, vectorL);
		Delta.calcDoubleDelta(delta, vectorL);

		return delta;
	}

	public static double centerFreq(int i, double samplingRate, double lowerFilterFreq, int numMelFilters) {
		double mel[] = new double[2];
		mel[0] = freqToMel(lowerFilterFreq);
		mel[1] = freqToMel(samplingRate / 2);

		// take inverse mel of:
		double temp = mel[0] + ((mel[1] - mel[0]) / (numMelFilters + 1)) * i;
		return inverseMel(temp);
	}

	public static double[] cepCoefficients(double f[], int numCepstra, int numFilters) {
		double cepc[] = new double[numCepstra];

		for (int i = 0; i < cepc.length; i++) {
			for (int j = 1; j <= numFilters; j++) {
				cepc[i] += f[j - 1] * Math.cos(Math.PI * i / numFilters * (j - 0.5));
			}
		}

		return cepc;
	}

	public static BufferedImage createImage(JPanel panel, int w, int h) {

		// int w = panel.getWidth();
		// int h = panel.getHeight();
		BufferedImage bi = new BufferedImage(w, h, BufferedImage.TYPE_INT_RGB);
		Graphics2D g = bi.createGraphics();
		panel.paint(g);
		return bi;
	}

	public static double[][] cutSpectrum(double[][] spectrum, float minFreq, float maxfreq, int fftWindowSize, int samplingRate) {
		int minFrequencyIndex = frequencyIndex(minFreq, fftWindowSize, samplingRate);
		int maxFrequencyIndex = frequencyIndex(maxfreq, fftWindowSize, samplingRate);

		double[][] cutSpectrum = new double[spectrum.length][maxFrequencyIndex - minFrequencyIndex + 1];

		for (int i = 0; i < spectrum.length; i++) {
			cutSpectrum[i] = Arrays.copyOfRange(spectrum[i], minFrequencyIndex, maxFrequencyIndex+1);
		}

		return cutSpectrum;
	}

	public static void exampleSet2Signal(double[] rebuiltSignal, ExampleSet es, Double fillerValueFormissingEntries) {

		MemoryExampleTable met = (MemoryExampleTable) es.getExampleTable();
		int numCol = met.getAttributeCount();
		int numRows = met.size();

		Attribute labelAtt = met.getAttribute(numCol - 1);

		for (int i = 0; i < numRows; i++) {
			int index = (int) met.getDataRow(i).get(labelAtt);
			String label = labelAtt.getMapping().mapIndex(index);
			int id = Integer.parseInt(label);
			Example e = es.getExample(i);
			// System.out.println(es.getExample(i)+"->"+signal[i]);
			for (Attribute a : e.getAttributes()) {
				Double value = e.getValue(a);
				if (value.equals(Double.NaN) && !fillerValueFormissingEntries.equals(Double.NaN))
					value = fillerValueFormissingEntries;

				rebuiltSignal[id] = value;
			}
		}
	}

	public static void exampleSet2Signal(double[] rebuiltSignal, ExampleSet es) {
		exampleSet2Signal(rebuiltSignal, es, null);
	}

	public static int[] fftBinIndices(double samplingRate, int frameSize, int numMelFilters, int numFequencies, float lowerFilterFreq, float upperFilterFreq) {
		int cbin[] = new int[numFequencies + 2];
		System.out.println("New Filter banks: " + numFequencies);
		cbin[0] = (int) Math.round(lowerFilterFreq / samplingRate * frameSize);
		cbin[cbin.length - 1] = frequencyIndex(upperFilterFreq, frameSize, (float) samplingRate);
		System.out.println("F0: " + lowerFilterFreq);
		for (int i = 1; i <= numFequencies; i++) {
			double fc = centerFreq(i, samplingRate, lowerFilterFreq, numMelFilters);
			System.out.println("F" + (i) + ": " + fc);
			cbin[i] = (int) Math.round(fc / samplingRate * frameSize);
		}

		System.out.println("F" + (cbin.length - 1) + ": " + upperFilterFreq);

		return cbin;
	}

	public static int[] fftBinIndices(double samplingRate, int frameSize, int numMelFilters, float lowerFilterFreq) {
		int cbin[] = new int[numMelFilters + 2];

		cbin[0] = (int) Math.round(lowerFilterFreq / samplingRate * frameSize);
		cbin[cbin.length - 1] = (int) (frameSize / 2);

		for (int i = 1; i <= numMelFilters; i++) {
			double fc = centerFreq(i, samplingRate, lowerFilterFreq, numMelFilters);

			cbin[i] = (int) Math.round(fc / samplingRate * frameSize);
		}

		return cbin;
	}

	public static double freqToMel(double freq) {
		return 2595 * log10(1 + freq / 700);
	}

	public static int frequencyIndex(float frequency, int fftSize, float samplingRate) {
		return Math.round(frequency * fftSize / samplingRate);
	}

	public static double[] generateSinSignal(int signalLength, float timeShift, float frequency) {
		// final float frequency = 0.3f;// 1f;

		double samples[] = new double[signalLength];
		float time = 0;
		for (int i = 0; i < samples.length; i++) {
			samples[i] = (float) Math.sin(2f * Math.PI * frequency * time);
			// time += 1f / (float) samplingRate;
			time += timeShift;
		}
		return samples;
	}

	public static double inverseMel(double x) {
		double temp = Math.pow(10, x / 2595) - 1;
		return 700 * (temp);
	}

	public static double log10(double value) {
		return Math.log(value) / Math.log(10);
	}

	public static double[] melFilter(double bin[], int cbin[], int numMelFilters) {
		double temp[] = new double[numMelFilters + 2];

		for (int k = 1; k <= numMelFilters; k++) {
			double num1 = 0, num2 = 0;

			for (int i = cbin[k - 1]; i <= cbin[k]; i++) {
				num1 += ((i - cbin[k - 1] + 1) / (cbin[k] - cbin[k - 1] + 1)) * bin[i];
			}

			for (int i = cbin[k] + 1; i <= cbin[k + 1]; i++) {
				num2 += (1 - ((i - cbin[k]) / (cbin[k + 1] - cbin[k] + 1))) * bin[i];
			}

			temp[k] = num1 + num2;
		}

		double fbank[] = new double[numMelFilters];
		for (int i = 0; i < numMelFilters; i++) {
			fbank[i] = temp[i + 1];
		}

		return fbank;
	}

	public static int recalculateMaxMelFilters(double samplingRate, int numMelFilters, float lowerFilterFreq, float maxFilterFreq) {
		int bestIndex = 1;
		for (int i = 1; i <= numMelFilters; i++) {
			double fc = centerFreq(i, samplingRate, lowerFilterFreq, numMelFilters);
			System.out.println("fc " + fc);
			if (fc > maxFilterFreq) {
				bestIndex = i;
				break;
			}
		}

		return bestIndex - 1;
	}

	
	public static double[] signalTimeLine(int signalLength, double samplingRate) {
		double time[] = new double[signalLength];
		Arrays.fill(time, Double.NaN);

		for (int i = 0; i < signalLength; i++) {
			time[i] = (double) i / (double) samplingRate;
		}
		System.out.println("time " + time[signalLength - 1] * samplingRate + " vs " + signalLength);
		return time;
	}

	public static float spectrumTime(float linearTime, float windowShiftTime) {
		return linearTime / windowShiftTime;
	}

	public static ExampleSet signal2ExampleSet(double[] signal) {
		BigSamplesTable samples = new BigSamplesTable();
		for (int k = 0; k < signal.length; k++) {
			samples.addSampleRow("" + k, signal[k]);
		}
		System.out.println("Example Set Created");
		return samples.generateExampleSet();
	}

	public static double[][] spectrogram(String name, double[] signal, int samplingRate, int windowshift, int frameslength, boolean display) throws Exception {
		SpectrogramCustom spec = new SpectrogramCustom(signal, samplingRate, Window.get(Window.HAMMING, frameslength), windowshift, frameslength, 640, 480);
		double[][] spectrum = spec.spectra.toArray(new double[spec.spectra.size()][]);
		if (display) {
			spec.showInJFrame(name, true, true);
			/*
			 * save spectrograms to files BufferedImage image = createImage(spec); ImageIO.write(ImageTools.toBufferedImage(image), "png", new File(name+".png"));
			 */
			// Thread.sleep(2000);
			// createImage(spec);
		}
		return spectrum;
	}

	public static void displaySpectrogram(double[][] spectrum, double[] signal, String name, int samplingRate, int windowshift, int frameslength) throws Exception {

		SpectrogramCustom spec = new SpectrogramCustom(signal, samplingRate, Window.get(Window.HAMMING, frameslength), windowshift, frameslength, 640, 480);

		spec.spectra = new ArrayList<double[]>();
		for (int i = 0; i < spectrum.length; i++) {
			spec.spectra.add(spectrum[i]);
		}
		spec.showInJFrame(name, true, true);

	}

	public static float spectrogramTimeFromIndex(int index, float windowShiftTime) {
		return index * windowShiftTime;
	}

	public static int spectrogramIndex(float linearTime, float windowShiftTime) {
		return (int) (linearTime / windowShiftTime);
	}

	public double[] averagepower;

	
	public int startStableTractIdx = -1;
	public int endStableTractIdx = -1;

	public double[] takeLongestStableTract(double[] signal, double valuedifftoleranceperc) {
		ArrayList<int[]> pairs = new ArrayList<int[]>();
		int idx1 = -1;

		int[] pair = null;
		// analyze the signal
		for (int i = 1; i < signal.length; i++) {
			// if there is not current range create it
			if (idx1 == -1) {
				idx1 = 1;
				pair = new int[2];
				pair[0] = i - 1;
				pair[1] = i - 1;
			}
			// if the current sample is similar to the previous, enlarge the range
			if (Math.abs(signal[i] - signal[i - 1]) / Math.max(signal[i], signal[i - 1]) <= valuedifftoleranceperc)
				pair[1] = i;
			// otherwise add the couple and reset
			else {
				idx1 = -1;
				pairs.add(pair);
			}
		}
		// if the last couple was reset, add the last interval
		if (idx1 > -1)
			pairs.add(pair);

		// find the longest pair
		int best = 0;
		int maxsize = 0;
		int k = 0;
		for (int[] setcouple : pairs) {
			int diff = setcouple[1] - setcouple[0];
			if (diff > maxsize) {
				maxsize = diff;
				best = k;
			}
			k++;
		}

		// take the longest range
		if (pairs.size() == 0) {
			pairs.add(new int[] { 0, 1 });
		}

		int[] bestcouple = pairs.get(best);
		// take the related slice of signal
		if (bestcouple[1]==bestcouple[0])
			bestcouple[1]=bestcouple[0]+1;
		double[] subsignal = new double[bestcouple[1] - bestcouple[0]];
		System.out.println("Longest range: from " + bestcouple[0] + " to " + bestcouple[1]);
		startStableTractIdx = bestcouple[0];
		endStableTractIdx = bestcouple[1];

		int l = 0;
		for (int i = bestcouple[0]; i < bestcouple[1]; i++) {
			subsignal[l] = signal[i];
			l++;
		}

		return subsignal;
	}



	// ###########Signal processing
	public static short[] silence(double durationInSec, double samplingFreq) {

		int nsamp = timeToSamples(durationInSec, samplingFreq);
		short[] d = new short[nsamp];
		for (int i = 0; i < nsamp; i++) {
			d[i] = 0;
		}
		return d;
	}

	public static double samplesToTime(int samples, double fs) {
		return (double) samples / fs;
	}

	public static int timeToSamples(double time, double fs) {
		return (int) Math.round(fs * time);
	}

	public static double[] featureTimesInSec(double windowShift, File audio) throws Exception {
		AudioBits bits = new AudioBits(audio);
		short[] signal = bits.getShortVectorAudio();
		bits.ais.close();
		float sfrequency = bits.getAudioFormat().getSampleRate();
		double total_audio_sec = (double) signal.length / (double) sfrequency;
		double t = 0;
		int ntimes = (int) Math.floor(total_audio_sec / windowShift);
		double times[] = new double[ntimes];
		int counter = 0;
		while (t < total_audio_sec) {
			if (counter < ntimes)
				times[counter] = t;
			t = t + windowShift;
			counter++;
		}
		return times;

	}

	public static double[] featureTimesInSec(double windowShift, int sfrequency, int signalLength) throws Exception {
		double total_audio_sec = (double) signalLength / (double) sfrequency;
		double t = 0;
		int ntimes = (int) Math.floor(total_audio_sec / windowShift);
		double times[] = new double[ntimes];
		int counter = 0;
		while (t < total_audio_sec) {
			if (counter < ntimes)
				times[counter] = t;
			t = t + windowShift;
			counter++;
		}
		return times;

	}
}
