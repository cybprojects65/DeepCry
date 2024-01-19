package it.cnr.speech.audiofeatures;

import java.io.File;

import it.cnr.workflow.utils.UtilsMath;

public class Energy {

	AudioBits bits;
	short[] signal;
	public double energyThr = 0.001;
	public float windowIns = 0.3f;
	public static String ENERGYSEGMENTATIONFILE = "energy_segmentation.csv";
	public static String DERIVATIVEFILE = "energy_derivative.csv";

	public float getWindowIns() {
		return windowIns;
	}

	public void setWindowIns(float windowIns) {
		this.windowIns = windowIns;
	}

	public double segmentSignal(File audioFile, File outputFolder) throws Exception {

		double[] normalisedEnergyCurve = energyCurve(windowIns, audioFile, false);
		double[] derivative = UtilsMath.derivative(normalisedEnergyCurve);
		double meanEnergy = UtilsMath.mean(normalisedEnergyCurve);
		// delete all segments
		if (!outputFolder.exists())
			outputFolder.mkdir();
				

		//Utils.writeSignal2File(derivative, new File(outputFolder,DERIVATIVEFILE));
		

		float sfrequency = bits.getAudioFormat().getSampleRate();
		int waveCounter = 0;
		int ntries = 0;
		int maxTries = 100;
		int minNumberOfWavesToFind = 3;
		double maxSNR = 0;
		
		while (ntries < maxTries) {
			maxSNR = 0;
			double time0 = 0;
			waveCounter = 0;
			//BufferedWriter bw = new BufferedWriter(new FileWriter(new File(outputFolder, ENERGYSEGMENTATIONFILE), true));

			for (int i = 1; i < normalisedEnergyCurve.length; i++) {
				if (derivative[i - 1] < 0 && normalisedEnergyCurve[i] < energyThr) {
					// normalisedEnergyCurve[i]=2;
					double time1 = ((double) (i) * Math.round(windowIns * sfrequency)) / (double) sfrequency;
					// System.out.println("Segment " + time1 + " : " +
					// normalisedEnergyCurve[i]);
					// save signal segment t0 t1
					int i0 = (int) (time0 * sfrequency);
					int i1 = (int) (time1 * sfrequency);
					short[] subsignal = new short[i1 - i0 + 1];
					for (int k = i0; k <= i1; k++) {
						subsignal[k - i0] = signal[k];
					}
					try {
						double SNR = 10*Math.log10(normalisedEnergyCurve[i-1]/normalisedEnergyCurve[i]);
						if (SNR>maxSNR)
							maxSNR=SNR;
						//System.out.println("SNR = "+SNR);
						File outputWaveFile = new File(outputFolder, "audio_segment_" + waveCounter + ".wav");
						AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(subsignal, outputWaveFile, bits.getAudioFormat());
						//bw.write(time0 + " " + time1 + " " + outputWaveFile.getName() + "\n");
					} catch (Exception e) {
						e.printStackTrace();
					}

					time0 = time1;
					waveCounter++;
				}
				if (derivative[i - 1] < 0) {
					//System.out.println("Energy with neg der " + " : " + normalisedEnergyCurve[i] + " vs " + energyThr);
				}
			}
			//bw.close();
			if (waveCounter > minNumberOfWavesToFind)
				ntries = maxTries;
			else {
				ntries++;
				System.out.println("Too few segments using energy threshold: " + energyThr + " mean energy " + meanEnergy);
				energyThr = energyThr + energyThr;
				System.out.println("Retrying segmentation using energy: " + energyThr);
			}

		} // end while

		if (waveCounter == 0)
			throw new Exception("Audio is too low for segmentation");
		
		System.out.println("SNR:"+maxSNR);
		
		return maxSNR;
		
	}

	public double[] energyCurve(float windowInMs, File audioFile) {
		return energyCurve(windowInMs, audioFile, true);
	}

	public double[] energyCurve(float windowInMs, File audioFile, boolean normalize) {
		// extract features
		bits = new AudioBits(audioFile);
		signal = bits.getShortVectorAudio();
		float sfrequency = bits.getAudioFormat().getSampleRate();
		double[] nrg = energyCurve(windowInMs, signal, sfrequency, normalize);
		// bits.deallocateAudio();
		return nrg;
	}

	public double[] energyCurve(float windowIns, short[] signal, float sfrequency, boolean normalize) {

		// initial energy
		int windowSamples = Math.round(windowIns * sfrequency);
		int steps = signal.length / windowSamples;

		// trace energy curve
		double[] energySignal = new double[steps];
		double maxEnergy = -Double.MAX_VALUE;
		for (int i = 0; i < steps; i++) {
			int currentIdx = i * windowSamples;
			short[] signalSlice = new short[windowSamples];
			for (int j = 0; j < windowSamples; j++) {
				signalSlice[j] = signal[currentIdx + j];
			}

			energySignal[i] = energy(signalSlice);

			if (energySignal[i] > maxEnergy)
				maxEnergy = energySignal[i];
		}

		if (normalize) {
			for (int i = 0; i < steps; i++) {
				energySignal[i] = energySignal[i] / maxEnergy;
			}
		}
		return energySignal;
	}
	
	public static double energy(short[] signal) {

		double energy = 0;
		for (int g = 0; g < signal.length; g++) {
			energy += signal[g] * signal[g];
		}
		energy = energy / (double) signal.length;
		return energy;
	}

}
