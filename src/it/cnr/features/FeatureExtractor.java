package it.cnr.features;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.io.FileUtils;

import it.cnr.speech.audiofeatures.Energy;
import it.cnr.speech.audiofeatures.PitchExtractor;
import it.cnr.speech.audiofeatures.SyllabicEnergyPitchExtractor;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.utilities.UtilsVectorMatrix;

public class FeatureExtractor {

	public double SNR;
	public WorkflowConfiguration config;
	
	public FeatureExtractor() {}
	
	public FeatureExtractor(WorkflowConfiguration config) {
		this.config = config;
	}
	
	public double getSNR() {
		return SNR;
	}

	public double[][] extractFeatureMatrix(File audioFile, boolean standardise){
		
		//extract energy
		List<double[]> featureList = new ArrayList<double[]>();
		double [] energyCurve = getEnergyFeatures(audioFile);
		System.out.println("Adding energy features");
		featureList.add(energyCurve);
		double [] pitchCurve = getPitchFeatures(audioFile);
		System.out.println("Adding pitch features");
		featureList.add(pitchCurve);
		
		//get min length of vectors
		int minElements = Integer.MAX_VALUE;
		for (int i=0;i<featureList.size();i++) {
			if (minElements>featureList.get(i).length)
				minElements=featureList.get(i).length;
		}
		System.out.println("Number of features: "+featureList.size());
		System.out.println("Number of samples per feature: "+minElements);
		System.out.println("Building feature matrix");
		//one feature vector for each row, one featue for each column
		int nrows = minElements;
		int ncols = featureList.size();
		double [][] featureMatrix =new double[nrows][ncols];
		for (int i=0;i<nrows;i++) {
			for (int j=0;j<ncols;j++) {
				featureMatrix[i][j] = featureList.get(j)[i];
			}
		}
		
		if (standardise) {
			System.out.println("Standardising the matrix...");
			UtilsVectorMatrix u = new UtilsVectorMatrix();
			featureMatrix = u.standardize(featureMatrix);
			
		}
		return featureMatrix;
	}
	
	public double [] getEnergyFeatures(File audioFile) {
		return new Energy().energyCurve(config.energyWindow4Analysis, audioFile, true);
	}
	
	public double [] getPitchFeatures(File audioFile) {
		PitchExtractor pitchExtr = new PitchExtractor();
		pitchExtr.setPitchWindowSec(config.pitchWindow4Analysis);
		pitchExtr.calculatePitch(audioFile.getAbsolutePath());
		Double [] pitchCurve = pitchExtr.pitchCurve;
		double pitch [] = new double[pitchCurve.length];
		for (int i=0;i<pitchCurve.length;i++) { 
			if (pitchCurve[i] == null || Double.isNaN(pitchCurve[i]) || Double.isInfinite(pitchCurve[i]))
				pitch [i]= 0;
			else
				pitch [i]= pitchCurve[i].doubleValue();
		}
		return pitch;
	}
	
	public File separateFilesBasedOnEnergy(File audioFile, float maxSilence) throws Exception {

		File outputFolder = new File(audioFile.getParentFile(), audioFile.getName().replace(".wav","_processing"));
		long t0 = 0;
		long t1 = 0;

		if (outputFolder.exists() && outputFolder.listFiles().length > 0)
			FileUtils.deleteDirectory(outputFolder);

		System.out.println("Temporary folder " + outputFolder);
		System.out.println("Segmenting audio based on silence energy (" + maxSilence + "s)");
		t0 = System.currentTimeMillis();
		Energy energy = new Energy();
		energy.setWindowIns(maxSilence);
		SNR = energy.segmentSignal(audioFile, outputFolder);
		System.out.println("Estimated SNR="+SNR);
		t1 = System.currentTimeMillis();
		System.out.println("Finished in " + ((float) (t1 - t0) / 1000f) + "s");

		return outputFolder;
	}

	public File extractEnergyPitchTimeSeries(File audioFile, float energyWindowSec, float pitchWindowSec)
			throws Exception {

		File outputFolder = new File(audioFile.getParentFile(),
				audioFile.getName() + "_energysegmentation_e" + energyWindowSec + "_p" + pitchWindowSec);
		if (outputFolder.exists() && outputFolder.listFiles().length > 0)
			FileUtils.deleteDirectory(outputFolder);
		outputFolder.mkdir();

		System.out.println("Extract Pitch and Energy from segments, using syllabic-size windows (energy:"
				+ energyWindowSec + "s; pitch:" + pitchWindowSec + "s)");

		long t0 = System.currentTimeMillis();
		SyllabicEnergyPitchExtractor epe = new SyllabicEnergyPitchExtractor(energyWindowSec, pitchWindowSec);
		epe.extractEnergyPitchOfWaveFiles(audioFile, outputFolder);
		long t1 = System.currentTimeMillis();
		System.out.println("Finished in " + ((float) (t1 - t0) / 1000f) + "s");

		return outputFolder;
	}

	public File extractEnergyPitchOfWaveFiles(File fs, float energyWindowSec, float pitchWindowSec) throws Exception {
		File outputFolder = new File(fs.getParentFile(),
				fs.getName() + "_energysegmentation_e" + energyWindowSec + "_p" + pitchWindowSec);
		if (outputFolder.exists() && outputFolder.listFiles().length > 0)
			FileUtils.deleteDirectory(outputFolder);
		outputFolder.mkdir();
		
		// File[] allClusterSegments = outputFolder.listFiles();
		File energyFile = new File(outputFolder, fs.getName() + "_" + SyllabicEnergyPitchExtractor.ENERGYFILE);
		if (energyFile.exists())
			energyFile.delete();

		File pitchFile = new File(outputFolder, fs.getName() + "_" + SyllabicEnergyPitchExtractor.PITCHFILE);
		if (pitchFile.exists())
			pitchFile.delete();

		BufferedWriter fw = new BufferedWriter(new FileWriter(energyFile, true));
		BufferedWriter fwP = new BufferedWriter(new FileWriter(pitchFile, true));

		if (fs.getName().endsWith(".wav") && fs.lastModified() > 0) {

			Energy nrg = new Energy();
			// EnergySegmenter nrg = new EnergySegmenter();
			// System.out.println("Extracting energy for " + fs.getName());
			double[] normalisedEnergyCurve = null;

			try {
				normalisedEnergyCurve = nrg.energyCurve(energyWindowSec, fs, true);
			} catch (Exception e) {
				System.out.println("Cannot extract energy from " + fs.getName() + " (insufficient data  < "
						+ energyWindowSec + " s )");
				fs.setLastModified(0);
			}

			for (int j = 0; j < normalisedEnergyCurve.length; j++) {
				double e = normalisedEnergyCurve[j];
				fw.write("" + e);
				if (j < normalisedEnergyCurve.length - 1)
					fw.write(",");
				else
					fw.write("\n");
			}
			if (normalisedEnergyCurve.length == 0)
				fw.write("\n");

			
			PitchExtractor syllabPitchExtr = new PitchExtractor();
			syllabPitchExtr.setPitchWindowSec(pitchWindowSec);
			
			double minPitch = 50;
			double maxPitch = 500;
			
			syllabPitchExtr.calculatePitch(fs.getAbsolutePath(),minPitch,maxPitch);

			Double[] pitch = syllabPitchExtr.pitchCurve;
			for (int j = 0; j < pitch.length; j++) {
				if (pitch[j] == null || Double.isNaN(pitch[j]) || Double.isInfinite(pitch[j]))
					fwP.write("0");
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

		fw.close();
		fwP.close();
		return outputFolder;
	}

	public File[] extractEnergyPitchTimeSeriesForAllWaves(File audioFolder, float energyWindowSec, float pitchWindowSec)
			throws Exception {
		long t0 = System.currentTimeMillis();

		List<File> allFolders = new ArrayList<>();
		File[] allFiles = audioFolder.listFiles();
		for (File audioFile : allFiles) {
			
			if (audioFile.getName().endsWith(".wav") && (audioFile.lastModified() != 0)) {
				File outputFolder = extractEnergyPitchOfWaveFiles(audioFile, energyWindowSec, pitchWindowSec);
				allFolders.add(outputFolder);
			}
		}

		long t1 = System.currentTimeMillis();
		System.out.println("All features extracted in " + ((float) (t1 - t0) / 1000f) + "s");

		File[] collectedFolders = new File[allFolders.size()];
		collectedFolders = allFolders.toArray(collectedFolders);
		return collectedFolders;
	}
	
	public int nFeaturesPerTime = 0;
	
	public double[][] chunkizeTimeSeries(File featurefolder,float energyWindow4Analysis, float pitchWindow4Analysis, float vectorframeinsec, float vectorframetimeshift) throws Exception{
		File [] featurefiles = featurefolder.listFiles();
		double [] energy = null;
		double [] pitch = null;
		nFeaturesPerTime = 2;
		
		for (File features:featurefiles) {
			
			if (features.getName().contains("energy") && features.getName().endsWith(".csv")) {
				String f = new String(Files.readAllBytes(features.toPath()));
				String [] fs = f.split(",");
				energy = new double[fs.length];
				for (int i = 0;i<fs.length;i++) {
					energy[i] = Double.parseDouble(fs[i]);
				}
				
			}else if (features.getName().contains("pitch") && features.getName().endsWith(".csv")){
				
				String f = new String(Files.readAllBytes(features.toPath()));
				String [] fs = f.split(",");
				pitch = new double[fs.length];
				for (int i = 0;i<fs.length;i++) {
					pitch[i] = Double.parseDouble(fs[i]);
				}
				
			}
			
			
		}
		
		float currenttime = 0;
		boolean collect = true;
		List<double[]> listOfFeatures = new ArrayList<>();
		int ncolumns = 0;
		
		while (collect) {
			
			float t0 = currenttime;
			float t1 = currenttime+vectorframeinsec;
			
			int energyIdx0 = (int)Math.floor(t0/energyWindow4Analysis); 
			int energyIdx1 = (int)Math.ceil(t1/energyWindow4Analysis)-1;
			
			int pitchIdx0 = (int)Math.floor(t0/pitchWindow4Analysis); 
			int pitchIdx1 = (int)Math.ceil(t1/pitchWindow4Analysis)-1;
			
			if (energyIdx1>(energy.length-1) || pitchIdx1>(pitch.length-1)) {
				collect = false;
				continue;
			}
			if (currenttime == 0) 
				ncolumns = (energyIdx1-energyIdx0+1)+(pitchIdx1-pitchIdx0+1);
			
			
			double [] vector = new double [ncolumns];
			
			for (int i=0;i<vector.length;i++) {
				if (i%2==0) {
					vector[i] = energy[(i/2)+energyIdx0];
				}else {
					vector[i] = pitch[((i-1)/2)+pitchIdx0];
				}
			}

			/*
			for (int i=0;i<vector.length;i++) {
				
				if (i<=(energyIdx1-energyIdx0))
					vector[i] = energy[i+energyIdx0];
				else
					vector[i] = pitch[(i-(energyIdx1-energyIdx0+1))+pitchIdx0];
			}
			*/
			listOfFeatures.add(vector);
			
			currenttime = currenttime+vectorframetimeshift;
		}
		
		double[][] featurematrix = new double[listOfFeatures.size()][ncolumns];
		featurematrix = listOfFeatures.toArray(featurematrix);
		
		System.out.println("Extracted "+featurematrix.length+" vectors of "+vectorframeinsec+"s duration from "+featurefolder.getName());
		
		return featurematrix;
		
	}
	
	
	public double[][] featureToTimeSeries(double  [] features, int vectorLength) throws Exception{
		int nrows = (int) (features.length/2);
		int ncols = vectorLength;
		double[][] timeseries = new double[nrows][ncols];
		
		for (int i = 0;i<nrows;i++) {
			
			for (int j=0;j<ncols;j++) {
				timeseries[i][j] = features[i*ncols+j];
			}
			
		}
		
		return timeseries;
		
	}
	
}
