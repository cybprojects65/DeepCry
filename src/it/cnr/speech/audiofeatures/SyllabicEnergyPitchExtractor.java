package it.cnr.speech.audiofeatures;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;

public class SyllabicEnergyPitchExtractor {

	public static String ENERGYFILE = "energy.csv";
	public static String PITCHFILE = "pitch.csv";
	public float energyWindowSec;
	public float pitchWindowSec;
	
	public SyllabicEnergyPitchExtractor(float energyWindowSec,float pitchWindowSec){
		this.energyWindowSec = energyWindowSec;
		this.pitchWindowSec = pitchWindowSec;
	}
	
	public void extractEnergyPitchOfWaveFiles(File folder) throws Exception {
		File[] allClusterSegments = folder.listFiles();
		File energyFile = new File(folder, ENERGYFILE);
		if (energyFile.exists())
			energyFile.delete();

		File pitchFile = new File(folder, PITCHFILE);
		if (pitchFile.exists())
			pitchFile.delete();

		BufferedWriter fw = new BufferedWriter(new FileWriter(energyFile, true));
		fw.write("filename,e0\n");
		BufferedWriter fwP = new BufferedWriter(new FileWriter(pitchFile, true));
		fwP.write("filename,is_question,f0\n");
		for (File fs : allClusterSegments) {
			if (fs.getName().endsWith(".wav") && fs.lastModified()>0) {

				Energy nrg = new Energy();
				//EnergySegmenter nrg = new EnergySegmenter();
				//System.out.println("Extracting energy for " + fs.getName());
				double[] normalisedEnergyCurve = null;
				
				try {
					normalisedEnergyCurve = nrg.energyCurve(energyWindowSec, fs, true);
				}catch(Exception e) {
					System.out.println("Cannot extract energy from "+fs.getName()+" (insufficient data  < "+energyWindowSec+" s ) - deleting");
					//e.printStackTrace();
					//System.out.println("Deleting file "+fs.getName()+" for insufficient length");
					fs.delete();
					fs.deleteOnExit();
					fs.setLastModified(0);
					continue;
				}
				
				fw.write(fs.getName() + ",");
				// System.out.println("N samples
				// "+normalisedEnergyCurve.length);
				for (int j = 0; j < normalisedEnergyCurve.length; j++) {
					double e = normalisedEnergyCurve[j];
					fw.write("" + e);
					if (j < normalisedEnergyCurve.length - 1)
						fw.write(",");
					else
						fw.write("\n");
				}
				if (normalisedEnergyCurve.length==0)
					fw.write("\n");
				
				PitchExtractor pitchExtr = new PitchExtractor();
				pitchExtr.calculatePitch(fs.getAbsolutePath());

				PitchExtractor syllabPitchExtr = new PitchExtractor();
				syllabPitchExtr.setPitchWindowSec(pitchWindowSec);
				boolean isq = syllabPitchExtr.isQuestionRough(fs.getAbsolutePath());
				//System.out.println("File " + fs.getName() + " is question? " + isq);
				fwP.write(fs.getName() + "," + isq + ",");

				Double[] pitch = pitchExtr.pitchCurve;
				for (int j = 0; j < pitch.length; j++) {
					if (pitch[j] == null || Double.isNaN(pitch[j]) || Double.isInfinite(pitch[j]))
						fwP.write(" ");
					else {
						double e = pitch[j];
						fwP.write("" + e);
					}
					if (j < pitch.length - 1)
						fwP.write(",");
					else
						fwP.write("\n");
				}
			}
		}

		fw.close();
		fwP.close();

	}
	
	
	public void extractEnergyPitchOfWaveFiles(File fs, File outputFolder) throws Exception {
		//File[] allClusterSegments = outputFolder.listFiles();
		File energyFile = new File(outputFolder, fs.getName()+"_"+ENERGYFILE);
		if (energyFile.exists())
			energyFile.delete();

		File pitchFile = new File(outputFolder, fs.getName()+"_"+PITCHFILE);
		if (pitchFile.exists())
			pitchFile.delete();

		BufferedWriter fw = new BufferedWriter(new FileWriter(energyFile, true));
		fw.write("filename,e0\n");
		BufferedWriter fwP = new BufferedWriter(new FileWriter(pitchFile, true));
		fwP.write("filename,is_question,f0\n");
		//for (File fs : allClusterSegments) 
		{
			if (fs.getName().endsWith(".wav") && fs.lastModified()>0) {

				Energy nrg = new Energy();
				//EnergySegmenter nrg = new EnergySegmenter();
				//System.out.println("Extracting energy for " + fs.getName());
				double[] normalisedEnergyCurve = null;
				
				try {
					normalisedEnergyCurve = nrg.energyCurve(energyWindowSec, fs, true);
				}catch(Exception e) {
					System.out.println("Cannot extract energy from "+fs.getName()+" (insufficient data  < "+energyWindowSec+" s ) - deleting");
					//e.printStackTrace();
					//System.out.println("Deleting file "+fs.getName()+" for insufficient length");
					fs.delete();
					fs.deleteOnExit();
					fs.setLastModified(0);
				}
				
				fw.write(fs.getName() + ",");
				// System.out.println("N samples
				// "+normalisedEnergyCurve.length);
				for (int j = 0; j < normalisedEnergyCurve.length; j++) {
					double e = normalisedEnergyCurve[j];
					fw.write("" + e);
					if (j < normalisedEnergyCurve.length - 1)
						fw.write(",");
					else
						fw.write("\n");
				}
				if (normalisedEnergyCurve.length==0)
					fw.write("\n");
				
				PitchExtractor pitchExtr = new PitchExtractor();
				pitchExtr.calculatePitch(fs.getAbsolutePath());

				PitchExtractor syllabPitchExtr = new PitchExtractor();
				syllabPitchExtr.setPitchWindowSec(pitchWindowSec);
				boolean isq = syllabPitchExtr.isQuestionRough(fs.getAbsolutePath());
				//System.out.println("File " + fs.getName() + " is question? " + isq);
				fwP.write(fs.getName() + "," + isq + ",");

				Double[] pitch = pitchExtr.pitchCurve;
				for (int j = 0; j < pitch.length; j++) {
					if (pitch[j] == null || Double.isNaN(pitch[j]) || Double.isInfinite(pitch[j]))
						fwP.write(" ");
					else {
						double e = pitch[j];
						fwP.write("" + e);
					}
					if (j < pitch.length - 1)
						fwP.write(",");
					else
						fwP.write("\n");
				}
			}
		}

		fw.close();
		fwP.close();

	}
	
	
	
}
