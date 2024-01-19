package it.cnr.features;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.workflow.utils.SignalProcessing;
import marytts.util.io.FileUtils;

public class CorpusCleaner {

	public static void deleteShortAudio(File audioFolder, float minimumLenghtinSec) throws Exception {

		File[] allFiles = audioFolder.listFiles();
		System.out.println("Deleting too short audio files");
		for (File audioFile : allFiles) {
			if (!audioFile.getName().endsWith(".wav"))
				continue;

			AudioBits bits = new AudioBits(audioFile);
			float sfrequency = bits.getAudioFormat().getSampleRate();
			int nsamples = bits.getShortVectorAudio().length;
			double audioseconds = SignalProcessing.samplesToTime(nsamples,sfrequency);
			bits.ais.close();

			if (audioseconds < minimumLenghtinSec) {

				System.out.println("Deleting audio " + audioFile.getName() + " Length " + audioseconds + "s < "
						+ minimumLenghtinSec + "s");
				audioFile.setLastModified(0);
				FileUtils.rename(audioFile.getAbsolutePath(), new File(audioFile.getAbsolutePath().replace(".wav", "_deleted.wav")).getAbsolutePath() );
				
				// FileUtils.forceDelete(audioFile);

			}else {
				System.out.println("Good audio " + audioFile.getName() + " Length " + audioseconds + "s > "
						+ minimumLenghtinSec + "s");
			}

		}

		System.out.println("Deleting too short audio files - Done");

	}

	public static File unifyToneUnitsIntoOneFile(File audioFolder, File audio) throws Exception {

		File[] allFiles = audioFolder.listFiles();
		System.out.println("Unifying tone units into one audio file");
		int totalSamples = 0;
		List<short[]> allGoodSignals = new ArrayList<>(); 
		for (File audioFile : allFiles) {
			if (audioFile.getName().endsWith("_deleted.wav"))
				continue;

			AudioBits bits = new AudioBits(audioFile);
			short [] signal = bits.getShortVectorAudio();
			bits.ais.close();
			totalSamples = totalSamples+signal.length;
			allGoodSignals.add(signal);
		}
		short[] unifiedSignal = new short[totalSamples];
		
		int idx = 0;
		for (int i = 0;i<allGoodSignals.size();i++) {
				
			short[] singlesignal = allGoodSignals.get(i);
			for (int j=0;j<singlesignal.length;j++) {
				unifiedSignal[idx] = singlesignal[j];
				idx = idx+1;
			}
		}
		
		File singleFile = new File(audioFolder,audio.getName().replace(".wav", "_shortened.wav"));
		
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(unifiedSignal, singleFile, new AudioBits(audio).getAudioFormat());
		
		System.out.println("Audio unified to "+singleFile.getAbsolutePath()+" - Done");
		
		return(singleFile);
	}
	
	public static void vector2LabFile(String [] annotations, double [] times, File labfile) throws Exception {

		BufferedWriter bw = new BufferedWriter(new FileWriter(labfile));

		int i = 0;
		for (String annotation: annotations) {
			double time0 = times[i];
			double time1 =  times[times.length-1];
			if (i<(times.length-1))
				time1 = times[i+1];
			String labLine = time0 + " " + time1 + " " + annotation;
			bw.write(labLine + "\n");
			i++;
		}

		bw.close();
	}
	
	public static void rawFormat2Lab(File rawfile, File labfile, double timeStep) throws Exception {

		String rawLine = new String(Files.readAllBytes(rawfile.toPath()));
		String[] values = rawLine.split(",");
		BufferedWriter bw = new BufferedWriter(new FileWriter(labfile));

		double time = 0;
		for (String val : values) {
			double time1 = Math.round((time + timeStep) * 100) / 100d;
			double valD = (double) Math.round(Double.parseDouble(val) * 100) / (double) 100;
			String labLine = time + " " + time1 + " " + valD;
			bw.write(labLine + "\n");
			time = time1;
		}

		bw.close();
	}

	public static void clusteringFormat2Lab(File rawfile, File labfile, double timeStep, double timeStepShift,
			List<Integer> clustersToReport) throws Exception {

		List<String> allLines = Files.readAllLines(rawfile.toPath());
		String[] rows = new String[allLines.size() - 1];
		int i = 0;
		for (String line : allLines) {
			if (i > 0) {
				String[] elements = line.split(",");
				int id = Integer.parseInt(elements[0]) - 1;
				double time0 = Math.round((id * timeStepShift) * 100) / 100d;
				double time1 = Math.round((time0 + timeStep) * 100) / 100d;// Math.round((time0+timeStepShift)*100)/100d;
				Integer cluster = Integer.parseInt(elements[2]);
				if (clustersToReport == null || clustersToReport.contains(cluster)) {
					String lineToAdd = time0 + " " + time1 + " " + cluster;
					rows[id] = lineToAdd;
				}
			}
			i++;
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(labfile));

		for (String val : rows) {
			if (val != null && val.length() > 0)
				bw.write(val + "\n");
		}

		bw.close();
	}

	public static void annotationVector2Lab(int[] annotations, File labfile, double timeStep, double timeStepShift,
			int exclude) throws Exception {

		String[] rows = new String[annotations.length];
		for (int i = 0; i < annotations.length; i++) {

			int id = annotations[i];
			if (exclude > -1 && exclude != id) {
				double time0 = Math.round((i * timeStepShift) * 100) / 100d;
				double time1 = Math.round((time0 + timeStep) * 100) / 100d;// Math.round((time0+timeStepShift)*100)/100d;

				String lineToAdd = time0 + " " + time1 + " " + id;
				rows[i] = lineToAdd;
			}
		}

		BufferedWriter bw = new BufferedWriter(new FileWriter(labfile));

		for (String val : rows) {
			if (val != null && val.length() > 0)
				bw.write(val + "\n");
		}

		bw.close();
	}

	public static void uniformLabels(File labelFile, File labeloutputFile) throws Exception {
		
		List<String> allLines = Files.readAllLines(labelFile.toPath());
		List<String> t0 = new ArrayList<>();
		List<String> t1 = new ArrayList<>();
		List<String> value = new ArrayList<>();
		String currentValue = null;
		List<String> newLines = new ArrayList<>();
		
		for (String line:allLines) {
			
			String [] elements = line.split(" ");
			String t00 = elements[0]; 
			String t10 = elements[1];
			String valueS = elements[2];
			
			if (currentValue!= null && !currentValue.equals(valueS)) {
				if (currentValue.length()>0) {
					String row = t0.get(0)+" "+t1.get(t1.size()-1)+" "+value.get(0);
					newLines.add(row);
				}
				t0 = new ArrayList<>();
				t1 = new ArrayList<>();
				value = new ArrayList<>();
				currentValue = valueS;	
			}else if (currentValue==null){
				currentValue = valueS;
			}
			
			t0.add(t00);
			t1.add(t10);
			value.add(valueS);
			
			
		}
		
		if (value.size()>0) {
			String row = t0.get(0)+" "+t1.get(t1.size()-1)+" "+value.get(0);
			newLines.add(row);
		}
		
		FileWriter fw = new FileWriter(labeloutputFile);
		for (int i=0;i<newLines.size();i++) {
			
			fw.write(newLines.get(i)+"\n");
			
		}
		
		fw.close();
				
		
	}
	
	public static void uniformLabelSegments(File labelFile, File labeloutputFile) throws Exception {
		
		List<String> allLines = Files.readAllLines(labelFile.toPath());
		
		double pt0 = -1;
		double pt1 = -1;
		String currentValue = null;
		
		List<String> newLines = new ArrayList<>();
		
		for (String line:allLines) {
			
			String [] elements = line.split(" ");
			double t00 = Double.parseDouble(elements[0]); 
			double t10 = Double.parseDouble(elements[1]);
			String valueS = elements[2];
			
			if (pt0<t00 && t00<pt1 && currentValue.equals(valueS)) {
				
				pt1 = t10;
				
			}else if (pt0==-1) {
				pt0 = t00;
				pt1 = t10;
				currentValue = valueS;
			}else {
				String row = pt0+" "+pt1+" "+currentValue;
				newLines.add(row);
				pt0 = t00;
				pt1 = t10;
				currentValue = valueS;
			}
			
		}
		
		if (pt0>-1) {
			String row = pt0+" "+pt1+" "+currentValue;
			newLines.add(row);
		}
		
		FileWriter fw = new FileWriter(labeloutputFile);
		for (int i=0;i<newLines.size();i++) {
			
			fw.write(newLines.get(i)+"\n");
			
		}
		
		fw.close();
				
		
	}

	public static void clusteringFormat2Lab(File rawfile, File labfile, double timeStep, double timeStepShift)
			throws Exception {
		clusteringFormat2Lab(rawfile, labfile, timeStep, timeStepShift, null);
	}

}
