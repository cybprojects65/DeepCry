package it.cnr.speech.audiofeatures;

import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.InputStream;
import java.nio.ByteBuffer;

import javax.sound.sampled.AudioFileFormat.Type;
import javax.sound.sampled.AudioFormat;
import javax.sound.sampled.AudioInputStream;
import javax.sound.sampled.AudioSystem;

public class AudioWaveGenerator {

	public static byte[] ShortToByte8(short[] buffer) {
		int N = buffer.length;
		float f[] = new float[N];
		float min = 0.0f;
		float max = 0.0f;
		for (int i = 0; i < N; i++) {
			f[i] = (float) (buffer[i]);
			if (f[i] > max)
				max = f[i];
			if (f[i] < min)
				min = f[i];
		}
		//float scaling = 1;
		float scaling = 1.0f + (max - min) / 256.0f; // +1 ensures we stay
														// within range and
														// guarantee no divide
														// by zero if sequence
														// is pure silence ...

		ByteBuffer byteBuf = ByteBuffer.allocate(N);
		for (int i = 0; i < N; i++) {
			byte b = (byte) (f[i] / scaling); /* convert to byte. */
			byteBuf.put(b);
		}
		return byteBuf.array();
	}

	public static void generateWave44khzFromSamples(short[] signal, File wavefile, AudioFormat format) throws Exception {
			//byte[] allBytes = ShortToByte(signal);
			ByteBuffer byteBuf = ByteBuffer.allocate(2*signal.length);
			
			for (short s : signal) byteBuf.putShort(s);
			byte[] allBytes = byteBuf.array();
			InputStream b_in = new ByteArrayInputStream(allBytes);
			AudioFormat nformat = //new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,format.getSampleRate(),16,format.getChannels(),format.getFrameSize(),format.getFrameRate(),format.isBigEndian());
					new AudioFormat(format.getSampleRate(), 16, 1, true, true);
			
			nformat = new AudioFormat(44100, 16, 1, true, true);
			
			AudioInputStream stream = new AudioInputStream(b_in, nformat,allBytes.length/format.getFrameSize());
			
			AudioSystem.write(stream, Type.WAVE, wavefile);
			b_in.close();
			stream.close();
	}
	
	public static void generateWaveFromSamplesWithSameFormat(short[] signal, File wavefile, AudioFormat format) throws Exception {
		//byte[] allBytes = ShortToByte(signal);
		ByteBuffer byteBuf = ByteBuffer.allocate(2*signal.length);
		
		for (short s : signal) byteBuf.putShort(s);
		byte[] allBytes = byteBuf.array();
		InputStream b_in = new ByteArrayInputStream(allBytes);
		AudioFormat nformat = //new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,format.getSampleRate(),16,format.getChannels(),format.getFrameSize(),format.getFrameRate(),format.isBigEndian());
				format; //new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		nformat = new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		AudioInputStream stream = new AudioInputStream(b_in, nformat,allBytes.length/format.getFrameSize());
		
		AudioSystem.write(stream, Type.WAVE, wavefile);
		b_in.close();
		stream.close();
	}
	
	public static void generateWaveFromSamplesWithSameFormat(double[] signal, File wavefile, AudioFormat format) throws Exception {
		short [] signalShort = AudioBits.doubleToShortVectorAudio(signal);
		generateWaveFromSamples(signalShort, wavefile, format);
	}
	
	public static void generateWaveFromSamples(short[] signal, File wavefile, AudioFormat format) throws Exception {
		//byte[] allBytes = ShortToByte(signal);
		ByteBuffer byteBuf = ByteBuffer.allocate(2*signal.length);
		
		for (short s : signal) byteBuf.putShort(s);
		byte[] allBytes = byteBuf.array();
		InputStream b_in = new ByteArrayInputStream(allBytes);
		AudioFormat nformat = //new AudioFormat(AudioFormat.Encoding.PCM_SIGNED,format.getSampleRate(),16,format.getChannels(),format.getFrameSize(),format.getFrameRate(),format.isBigEndian());
				new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		nformat = new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		AudioInputStream stream = new AudioInputStream(b_in, nformat,allBytes.length/format.getFrameSize());
		
		AudioSystem.write(stream, Type.WAVE, wavefile);
		b_in.close();
		stream.close();
}
	
	public static void generateWaveFromSamples(double[] signal, File wavefile, AudioFormat format) throws Exception {

		ByteBuffer byteBuf = ByteBuffer.allocate(8*signal.length);
		
		for (double s : signal) 
			byteBuf.putDouble(s);
		
		byte[] allBytes = byteBuf.array();
		InputStream b_in = new ByteArrayInputStream(allBytes);
		AudioFormat nformat = 
				new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		nformat = new AudioFormat(format.getSampleRate(), 16, 1, true, true);
		
		AudioInputStream stream = new AudioInputStream(b_in, nformat,allBytes.length/format.getFrameSize());
		
		AudioSystem.write(stream, Type.WAVE, wavefile);
		b_in.close();
		stream.close();
	}
	

}
