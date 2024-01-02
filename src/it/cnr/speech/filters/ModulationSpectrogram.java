package it.cnr.speech.filters;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.Arrays;

import javax.sound.sampled.AudioFormat;

import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;

public class ModulationSpectrogram {

	public double[][] modulationSpectrogram = null;
	public static double windowSize = 0.250; // in s
	public static double windowShift = 0.0125; // in s
	public static double cutoffFreq = 28;// Hz
	public static int reductionFactor = 100;
	public static boolean saveSteps = false;
	public static double maxdB = 30;
	public static double mindB = -30;
	
	public void calcMS(File audio, File output, boolean saturate, boolean addDeltas, int numberOfMSFeatures, double maxFrequency) throws Exception {

		SignalConverter sc = new SignalConverter();
		sc.getSignal(audio);

		/*
		 * LowPassFilterDynamic lptest = new LowPassFilterDynamic(cutoffFreq,
		 * sc.samplingRate, 0.01); //LowPassFilterStatic lowpass = new
		 * LowPassFilterStatic(cutoffFreq, sc.samplingRate); double [] lptestsig =
		 * lptest.applyFilter(sc.signal); File stepCheckFile1 = new
		 * File(audio.getAbsolutePath().replace(".wav", "_lowp28orig.wav"));
		 * AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(lptestsig,
		 * stepCheckFile1, new AudioBits(audio).getAudioFormat());
		 * 
		 * System.exit(0);
		 */
		
		for (int k = 0; k < numberOfMSFeatures; k++) {
			int melFilterCounter = k;
			// I will use the original Greenwood filterbank from "Critical bandwidth and the
			// frequency coordinates of the basilar membrane." simulated as Mel triangular
			// filters
			// using 15 Mel filters
			System.out.println("#RUNNING Mel filter n. " + melFilterCounter);

			GreenwoodFilterBank mel = new GreenwoodFilterBank(melFilterCounter, sc.samplingRate, 0.01,numberOfMSFeatures,maxFrequency);
			double[] mel_signal = mel.applyFilter(sc.signal);
			File stepCheckFile = null;

			if (saveSteps) {
				stepCheckFile = new File(audio.getAbsolutePath().replace(".wav", "_mel_" + melFilterCounter + ".wav"));
				AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(mel_signal, stepCheckFile,
						new AudioBits(audio).getAudioFormat());
			}
			System.out.println("Mel filter n. " + melFilterCounter + " done");

			System.out.println("Rectification");
			double[] hwr_signal = HalfWaveRectification.applyHalfWaveRectification(mel_signal);
			if (saveSteps) {
				stepCheckFile = new File(stepCheckFile.getAbsolutePath().replace(".wav", "_rect.wav"));
				AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(hwr_signal, stepCheckFile,
						new AudioBits(audio).getAudioFormat());
			}
			System.out.println("Rectification done");

			System.out.println("Low-pass of " + 28 + "Hz");
			LowPassFilterDynamic lowpass = new LowPassFilterDynamic(cutoffFreq, sc.samplingRate, 0.01);
			
			double[] lowpassed_signal = lowpass.applyFilter(hwr_signal);
			System.out.println("Low-pass OK");
			if (saveSteps) {
				stepCheckFile = new File(stepCheckFile.getAbsolutePath().replace(".wav", "_lowp.wav"));
				AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(lowpassed_signal, stepCheckFile,
						new AudioBits(audio).getAudioFormat());
			}

			System.out.println("Downsampling");
			double[] lowpassed_subsampled_signal = DownSampler.downsample(lowpassed_signal, reductionFactor);
			int samplingRateReduced = sc.samplingRate / reductionFactor;
			AudioFormat downaf = new AudioFormat(samplingRateReduced, 16, 1, true, true);
			if (saveSteps) {
				stepCheckFile = new File(stepCheckFile.getAbsolutePath().replace(".wav", "_downsp.wav"));
				AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(lowpassed_subsampled_signal, stepCheckFile,
						downaf);
			}
			System.out.println("Downsampling OK");

			System.out.println("Normalising");
			double[] lowpassed_subsampled_normalised_signal = SignalConverter
					.normaliseByAverageEnvelopeLevel(lowpassed_subsampled_signal);

			if (saveSteps) {
				stepCheckFile = new File(stepCheckFile.getAbsolutePath().replace(".wav", "_norm.wav"));
				AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(lowpassed_subsampled_normalised_signal,
						stepCheckFile, downaf);
			}
			System.out.println("Normalising OK");

			System.out.println("FFT");
			SignalConverter scReduced = new SignalConverter();
			double[][] spectrogram = scReduced.shortTermFFT(lowpassed_subsampled_normalised_signal, samplingRateReduced,
					windowSize, windowShift);
			System.out.println("Spectrogram FFT done");
			float minPossibleFreq = 4;
			float maxPossibleFreq = 4;

			System.out.println("Cut-off spectrum to 4Hz");
			//public static double[][] cutSpectrum(double[][] spectrum, float minFreq, float maxfreq, int fftWindowSize, int samplingRate) {
			
			double[][] spectrogram_4hz = SignalConverter.cutSpectrum(spectrogram, 
					minPossibleFreq, maxPossibleFreq,
					scReduced.windowSizeSamples, samplingRateReduced);
			double[] spectrogram_4hz_slice = new double[spectrogram_4hz.length];
			
			
			//report magnitude in Dbs and saturate dBs between 30 and -30 dBs
			
			for (int i = 0; i < spectrogram_4hz_slice.length; i++) {
					spectrogram_4hz_slice[i] = Math.pow(spectrogram_4hz[i][0], 2);
					
					if (saturate) {
					//convert to dB and saturate over 30 and under 0
						if (spectrogram_4hz_slice[i] == 0)
							spectrogram_4hz_slice[i] = -30;
						else {
							spectrogram_4hz_slice[i] = 10 * Math.log10(spectrogram_4hz_slice[i]);
							if (spectrogram_4hz_slice[i]<=0)
								spectrogram_4hz_slice[i] = -30;
							else if (spectrogram_4hz_slice[i]>30)
								spectrogram_4hz_slice[i] = 30;
						}
					} 
					
				}
			
			
			System.out.println("Cut-off spectrum OK");

			if (modulationSpectrogram == null)
				modulationSpectrogram = new double[numberOfMSFeatures][spectrogram_4hz_slice.length];

			modulationSpectrogram[numberOfMSFeatures-k-1] = spectrogram_4hz_slice;

			
		}

		//double [][] modulationSpectrogramD = addDelta(modulationSpectrogram);
		if (addDeltas)
			modulationSpectrogram = addDelta(modulationSpectrogram);
		if (output != null)
			save(output);

	}

	public static double[][] addDelta(double[][] features) throws Exception {
		int frames = features[0].length;
		double[][] delta = new double[features.length*2][frames];

		for (int i = 0;i<features.length;i++) {
			
			for (int j = 0;j<frames;j++) {
				delta[i+features.length][j] = features[i][j];
				double curr = features[i][j];
				double nex = curr;
				
				if (j< (features.length-1) ) {
					nex = features[i][j+1];
				}
				
				double D = nex-curr;
				delta[i][j] = D;
			}
			
		}
		
		
		return delta;
	}
	
	public void save(File outputFile) throws Exception {
		BufferedWriter bw = new BufferedWriter(new FileWriter(outputFile));

		for (int i = 0; i < modulationSpectrogram.length; i++) {

			String modspecSingle = Arrays.toString(modulationSpectrogram[i]);
			bw.write(modspecSingle.replace("[", "").replace("]", "").replace(" ", "") + "\n");

		}

		bw.close();

	}

}
