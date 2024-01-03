package it.cnr.clustering;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.cnr.features.CorpusCleaner;
import it.cnr.features.Utils;
import it.cnr.speech.filters.SignalConverter;
import it.cnr.workflow.configuration.Configuration;

public class ClassificationManager extends ClusteringManager{
	
	public ClassificationManager(Configuration config, File audio) {
		super(config, audio);
	}
	
	public void clusterFeatures(double [][] features, File outputFolder, int minClusters, int maxClusters, double entropyThreshold, double lowestEntropyThreshold) throws Exception{
		System.out.println("Starting multi-k-means clustering");
		MultiKMeans clusterer = new MultiKMeans();
		File clusterFile = clusterer.clusterFeatures(features, outputFolder, config.minNFeaturesInCluster,minClusters,maxClusters);
		clusteredFeatures = new HashMap<Integer, Cluster>();
		vectorID2ClusterID = new HashMap<Integer, Integer>();
		
		List<String> all = Files.readAllLines(clusterFile.toPath());
		int i = 0;
		for (String line : all) {
			if (i > 0) {

				String[] row = line.split(",");
				int id = Integer.parseInt(row[0]) - 1;
				Integer cluster = Integer.parseInt(row[2]);
				Boolean isOutlier = Boolean.parseBoolean(row[3]);
				double [] vector = features[id];
				Cluster c  = clusteredFeatures.get(cluster);
				if (c == null) {
					c = new Cluster(cluster, isOutlier);
				}
				
				c.add(vector, id);
				clusteredFeatures.put(cluster, c);
				vectorID2ClusterID.put( id, cluster);
			}

			i++;
		}
		System.out.println("Number of clusters identified: "+clusteredFeatures.keySet().size()) ;
		
		centroids = new HashMap<Integer, double[]>();
		centroid_interpretations = new HashMap<Integer, String>();
		System.out.println("Centroid interpretation");
		
		for (Integer id:clusteredFeatures.keySet()) {
			Cluster c = clusteredFeatures.get(id);
			System.out.println("Cluster "+id+" abundance: "+c.nElements) ;
			double [] centroid = c.calcCentroid();
			double entropy = Utils.roundDecimal(SignalConverter.calculateSpectralEntropy(centroid),2);
			double energy = Utils.roundDecimal(Utils.mean(centroid),2);
			double indicator =Utils.roundDecimal(energy*entropy,2);
			//String entropyInterpretation = ""+entropy;
			//
			String centroidIndicator = indicator+" ["+entropy+"/"+energy+"]"+" ("+id+")";
			
			System.out.println("Cluster "+id+" indicator: "+centroidIndicator) ;
			String entropyInterpretation = " ";
			if (indicator >lowestEntropyThreshold && indicator < entropyThreshold) {
			//if (energy < 17)
				entropyInterpretation = "Anomalous "+"("+id+")";
			}else {
				//entropyInterpretation = "N "+"("+id+")";
			}
			
			String centroidInterpretation = entropyInterpretation;
			
			centroid_interpretations.put(id, centroidInterpretation);
			centroids.put(id, centroid);
			System.out.println("Centroid with ID "+id+": interpreted as "+centroidInterpretation);
		}
		
	}
	
	public File toLabCTC(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec, double minLabelTimeSec) throws Exception{
		
		//toLab(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec, String [] labels) throws Exception{
		times = Utils.featureTimesInSec(windowLengthSec,samplingFrequency,signalLength);

		int nfeatures = vectorID2ClusterID.keySet().size();
		
		labels = new String[times.length];
		
		for (int i=0;i<nfeatures;i++) {
			Integer id = i;
			Integer clusterID = vectorID2ClusterID.get(id);
			String interpretation = " ";
			if (clusterID!=null) {
				interpretation = centroid_interpretations.get(clusterID);
			}
			labels[i] = interpretation; 
		}
		
		//merge consecutive labels
		boolean merge = true;
		if (merge) {
		int i=0;
		List<Double> timesNew = new ArrayList<Double>();
		List<String> labelsNew = new ArrayList<String>();
		String prevLab = "S";
		
		while(i<times.length) {
			String currentLab = labels[i];
			if (currentLab==null)
				currentLab=" ";
			if (!prevLab.equals(currentLab)) {
					timesNew.add(times[i]);
					labelsNew.add(currentLab);
					prevLab=currentLab;
			}else {
					
					
				}
			i++;
			}
		
		timesNew.add(times[times.length-1]);
		labelsNew.add(" ");
		
		i=0;
		int ntimes = timesNew.size();
		double lastTime = times[times.length-1];
		
		while(i<ntimes) {
			Double currTime = timesNew.get(i);
			Double nextTime = lastTime;
			if (i<(ntimes-1)) {
				nextTime = timesNew.get(i+1);
			}
			if ((nextTime-currTime)<minLabelTimeSec) {
				labelsNew.set(i, " ");
			}else {
					System.out.println("OK: SUITABLE RANGE FOUND: "+(nextTime-currTime)+"s");
			}
			
			i++;
		}
				
		times = new double[timesNew.size()];
		labels = new String[timesNew.size()];
        for (int j = 0; j < times.length; j++) {
            times[j] = timesNew.get(j).doubleValue();
            labels[j] = labelsNew.get(j);
        }
        
		}
		CorpusCleaner.vector2LabFile(labels, times,outputFile);
		return outputFile;
	}

}
