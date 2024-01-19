package it.cnr.speech.audiofeatures;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.IOException;
import java.io.InputStream;
import java.net.URL;
import java.util.Vector;

import javax.sound.sampled.AudioFileFormat;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.DataLine;
import javax.sound.sampled.TargetDataLine;

public class VAD {

	AudioInputStream inputStream;
	
	public File recordedfile;
	public boolean bargein = false;
	public float maxInitialSilence = 20; //seconds of silence without uttering a word
	public float maxSilenceAfterWord = -0.3f; //seconds of silence after a word
	public float maxSpeechTime = 10; //max speech in seconds
	public float sensitivity = 100;
	public TargetDataLine line;
	public int bytesToRead = 1600*5;//1600 = 100ms // old value 160 = 0.1s;
	public double EnergyThr = Math.pow(10, 8);
	public boolean processend = false;
	public AudioFormat audioFormat ;
	
	private Vector<byte[]> packets = new Vector<byte[]>();
	private int bargeIndex = 0;
	private int endIndex = -1;
	private boolean acquiring = true;
	
	public synchronized void changeAcquisitionStatus(boolean status) {
		acquiring = status;
	}

	public synchronized boolean getAcquisitionStatus() {
		return acquiring;
	}

	public void setMaxSilence(int packs) {
		maxSilenceAfterWord = packs;
	}

	public void init(File recordingFile) throws Exception {
		recordedfile = recordingFile;
		
		//InputStream in = ClassLoader.getSystemClassLoader().getResourceAsStream("format.wav");
		
		URL formatFile = ClassLoader.getSystemClassLoader().getResource("format.wav");
		//AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(new File(formatFile.toURI()));
		
		//AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(in);
		AudioInputStream audioInputStream = AudioSystem.getAudioInputStream(formatFile);
		audioFormat = audioInputStream.getFormat();
		System.out.println("Recording in " + audioFormat);
		DataLine.Info info = new DataLine.Info(TargetDataLine.class,
				audioFormat);
		line = (TargetDataLine) AudioSystem.getLine(info);
		line.open(audioFormat);
		inputStream = new AudioInputStream(line);
	}

	public boolean catchAudio(File recordingFile) throws Exception {
		init(recordingFile);
		System.out.println("Mic initialised..");
		PacketsWriter sp = new PacketsWriter();
		Thread t = new Thread(sp);
		t.start();
		System.out.println("Listener started. Analysing..");
		analyze();
		System.out.println("Audio Analysed. Waiting for flushing..");
		wait4ProcessEnd();
		System.out.println("Process finished.");
		return bargein;
	}

	public void wait4ProcessEnd() throws Exception {

		while (!processend) {
			Thread.sleep(200);
		}
	}

	public void closeaudio(){
    	System.out.println("drain");
		line.drain();
		System.out.println("close");
		line.close();
		System.out.println("stop");
		line.stop();
    }
	
	public void analyze() {

		try {
			double initialEnergy = -1;
			boolean lastbarged = false;
			double timecounter = 0;
			long time = 0;
			
			while ((packets == null) || (packets.size() == 0)) {
			}

			int frameSize = inputStream.getFormat().getFrameSize();
			int j = 0;

			while ((getAcquisitionStatus()) && (time < maxSpeechTime)) {
				while (packets.size() < j + 1) {
				}

				byte[] audiobytes = packets.get(j);
				int length = audiobytes.length / frameSize;

				double[] samples = new double[length];
				int db = 0;
				int start = 0;
				for (int g = 0; g < length * 2; g = g + frameSize) {
					db = (int) audiobytes[g + frameSize - 1];
					for (int b = frameSize - 2; b >= 0; b--)
						db = db << 8 | ((int) audiobytes[g + b] & 0xff);
					samples[start++] = ((double) db);
				}

				double energy = 0;
				for (int g = 0; g < length; g++) {
					energy += samples[g] * samples[g];
				}
				double timelength = (samples.length/audioFormat.getSampleRate());
				
				energy = energy / samples.length;
				if ((initialEnergy == -1) && (energy < EnergyThr)) {
					
					initialEnergy = energy;
					//System.out.println("N. samples of captured audio "+samples.length+" = "+timelength+"s");
					//System.out.println("Initial energy : "+initialEnergy);
				}
				
				
				
				//if ((energy > (initialEnergy + sensitivity * initialEnergy))
				if ((energy / initialEnergy) > sensitivity) 
				//		|| (energy > EnergyThr)) 
				{
					
					System.out.println("energy: "+(energy/initialEnergy));
					if (!bargein) {
						if (initialEnergy == -1)
							initialEnergy = EnergyThr;
						bargein = true;
						bargeIndex = j;
					}
					lastbarged = true;
				}
				else {
					if (lastbarged) {
						lastbarged = false;
						timecounter = 0;
					} else {
						timecounter=timecounter+timelength;
					}
					//System.out.println("Time counter "+timecounter);
					if ((bargein) && (timecounter > maxSilenceAfterWord)) {
						//System.out.println("Max silence after word reached");
						changeAcquisitionStatus(false);
						endIndex = j;
						break;
					} else if (timecounter > maxInitialSilence) {
						//System.out.println("Nothing was uttered");
						changeAcquisitionStatus(false);
						break;
					}
				}

				j++;
				time += timelength;
			}

			if (endIndex == -1)
				endIndex = j - 1;

			changeAcquisitionStatus(false);

		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public byte[] generateBytes() {

		byte[] audio = new byte[((endIndex - bargeIndex + 1) * bytesToRead)];
		int k = 0;
		try {
			for (int i = bargeIndex; i <= endIndex; i++) {
				byte[] audiopack = packets.get(i);
				for (int j = 0; j < bytesToRead; j++) {
					audio[k] = audiopack[j];
					k++;
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error in acquisition packet n. " + k);
		}
		return audio;
	}

	public void saveBytes(File filename, byte[] data) {
		try {
			ByteArrayInputStream bais = new ByteArrayInputStream(data);
			int frameSize = inputStream.getFormat().getFrameSize();
			int length = data.length / frameSize;
			AudioInputStream ais = new AudioInputStream(bais,
					inputStream.getFormat(), length);
			AudioSystem.write(ais, AudioFileFormat.Type.WAVE, recordedfile);
			processend = true;
		} catch (Exception e) {
			e.printStackTrace();
		}
	}


	public class PacketsWriter implements Runnable {

		public void run() {

			changeAcquisitionStatus(true);
			int bytesRead = 0;
			byte[] buffer = new byte[bytesToRead];

			line.start();
			int i = 0;
			while (getAcquisitionStatus() && (bytesRead >= 0)) {
				try {
					bytesRead = read(buffer, bytesToRead, inputStream);
					packets.add(i, buffer.clone());
					i++;
				} catch (IOException e) {
					System.out.println("Error in bytes acquisition");
				}

			}
			System.out.println("closing audio");
			closeaudio();
			System.out.println("saving bytes");
			saveBytes(recordedfile, generateBytes());
			buffer = null;
		}

		protected int read(byte buffer[], int bytesToRead, AudioInputStream ais)
				throws IOException {
			int nBytesRead = 0;
			int nBytesTotalRead = 0;

			while (nBytesRead != -1 && bytesToRead > 0) {
				nBytesRead = ais.read(buffer, nBytesTotalRead, bytesToRead);
				if (nBytesRead != -1) {
					bytesToRead -= nBytesRead;
					nBytesTotalRead += nBytesRead;
				}
			}
			return nBytesTotalRead;
		}
	}

}