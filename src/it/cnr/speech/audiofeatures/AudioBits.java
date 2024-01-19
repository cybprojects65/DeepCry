package it.cnr.speech.audiofeatures;

import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;
import javax.sound.sampled.UnsupportedAudioFileException;

import java.io.*;

/**
 * DEFAULT_SAMPLE_RATE = 8000; DEFAULT_FRAME_RATE = 8000; DEFAULT_CHANNELS = 1;
 * DEFAULT_FRAME_SIZE = 1; DEFAULT_SIZE_IN_BITS = 8; DEFAULT_ENCODING = "ALAW";
 * </PRE>
 */
public class AudioBits {

	private File file = null;
	private double samples[] = null;
	private short shortSamples[] = null;
	public AudioInputStream ais = null;

	public AudioBits(File f) {
		file = f;
		loadFile();
	}

	
	public AudioFormat getAudioFormat(){
		return ais.getFormat();
	}
	
	private void loadFile() {
		try {
			ais = AudioSystem.getAudioInputStream(file);
		}catch (Exception e) {
			e.printStackTrace();
			System.out.println("Error while reading audio stream from "
					+ file.getAbsolutePath());
		}
	}

	/**
	 * Generates a short numbers vector containing the audio samples
	 */
	public short[] getShortVectorAudio() {
		if (samples == null) {
			this.samples = getDoubleVectorAudio();
		}

		int length = this.samples.length;

		this.shortSamples = new short[length];
		for (int i = 0; i < length; i++) {
			Double tmp = this.samples[i];
			shortSamples[i] = tmp.shortValue();
		}
		return this.shortSamples;
	}
	
	public static short[] doubleToShortVectorAudio(double [] samples) {
		
		int length = samples.length;

		short [] shortSamples = new short[length];
		for (int i = 0; i < length; i++) {
			Double tmp = samples[i];
			shortSamples[i] = tmp.shortValue();
		}
		return shortSamples;
	}
	
	/**
	 * Generates a double numbers vector containing the audio samples
	 * 
	 */
	public double[] getDoubleVectorAudio() {
		int bytesRead = 0;
		AudioFormat source = ais.getFormat();
		long len = (long) ais.getFrameLength();
		int frameSize = source.getFrameSize();
		int bytesToRead = (int) len * frameSize;
		byte[] buffer = new byte[bytesToRead];
		try {
			bytesRead = read(buffer, bytesToRead);
		} catch (IOException e) {
			e.printStackTrace();
			System.out
					.println("Error during audio reading");
		}

		int window = 64; 
		int parts = (bytesRead / frameSize) / window;

		int length = (parts * window);
		samples = new double[length];

		int db = 0;
		int start = 0;

		// the bytes of one sample value are in little-Endian order
		for (int i = 0; i < length * 2; i = i + frameSize)
		{
			// first byte will be converted to int with respect to the sign
			db = (int) buffer[i + frameSize - 1];

			// combine the bytes of the sample (unsigned bytes) in reversed
			// order
			for (int b = frameSize - 2; b >= 0; b--)
				db = db << 8 | ((int) buffer[i + b] & 0xff);

			// convert to double value 
			samples[start++] = ((double) db);
		}
		
		buffer = null;
		return samples;
	}
	
	
	private int read(byte buffer[], int bytesToRead) throws IOException {
		int nBytesRead = 0;
		int nBytesTotalRead = 0;

		while (nBytesRead != -1 && bytesToRead > 0)
		{
			nBytesRead = ais.read(buffer, nBytesTotalRead, bytesToRead);
			if (nBytesRead != -1) {
				bytesToRead -= nBytesRead;
				nBytesTotalRead += nBytesRead;
			}
		}
		return nBytesTotalRead;
	}

	
	public void deallocateAudio() {
		try {
			ais.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		this.ais = null;
	}


	public static byte [] getAudioDataBytes(byte [] sourceBytes, AudioFormat audioFormat) throws UnsupportedAudioFileException, IllegalArgumentException, Exception {
	    if(sourceBytes == null || sourceBytes.length == 0 || audioFormat == null){
	        throw new IllegalArgumentException("Illegal Argument passed to this method");
	    }

	    try (final ByteArrayInputStream bais = new ByteArrayInputStream(sourceBytes);
	         final AudioInputStream sourceAIS = AudioSystem.getAudioInputStream(bais)) {
	        AudioFormat sourceFormat = sourceAIS.getFormat();
	        AudioFormat convertFormat = new AudioFormat(AudioFormat.Encoding.PCM_SIGNED, sourceFormat.getSampleRate(), 16, sourceFormat.getChannels(), sourceFormat.getChannels()*2, sourceFormat.getSampleRate(), false);
	        try (final AudioInputStream convert1AIS = AudioSystem.getAudioInputStream(convertFormat, sourceAIS);
	             final AudioInputStream convert2AIS = AudioSystem.getAudioInputStream(audioFormat, convert1AIS);
	             final ByteArrayOutputStream baos = new ByteArrayOutputStream()) {
	            byte [] buffer = new byte[8192];
	            while(true){
	                int readCount = convert2AIS.read(buffer, 0, buffer.length);
	                if(readCount == -1){
	                    break;
	                }
	                baos.write(buffer, 0, readCount);
	            }
	            return baos.toByteArray();
	        }
	    }
	}
	
}
