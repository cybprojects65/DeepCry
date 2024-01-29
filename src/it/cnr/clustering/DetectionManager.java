package it.cnr.clustering;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.cnr.features.CorpusCleaner;
import it.cnr.workflow.utils.SignalProcessing;

public class DetectionManager {
	
	public HashMap<Integer, Cluster> clusteredFeatures; 
	public HashMap<Integer, Integer> vectorID2ClusterID; 
	public HashMap<Integer, double [] > centroids;
	public HashMap<Integer, String> centroid_interpretations;
	//public WorkflowConfiguration config;
	public File audio;
	public static String HIGH_EP = "H";
	public boolean highriskclusterfound = false;
	public String labels [];
	public double times[];
	
	public DetectionManager(File audio) {
		this.audio = audio;
	}
	
	public void detectHighValuedFeatures(double [][] features, File outputFolder, int minClusters, int maxClusters, double[] thresholdPerFeature) throws Exception{
		System.out.println("Starting multi-k-means clustering");
		
		//XMeansClustering clusterer = new XMeansClustering();
		//File clusterFile = clusterer.clusterFeatures(features, outputFolder, config.minNFeaturesInCluster,1,10);
		//MultiKMeans clusterer = new MultiKMeans();
		//File clusterFile = clusterer.clusterFeatures(features, outputFolder, config.minNFeaturesInCluster,minClusters,maxClusters);
		
		QuickKMeans cl = new QuickKMeans();
		File clusterFile = cl.kMeans(features, maxClusters , outputFolder);
		
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
		//double [] means = Utils.columnMeans(features);
		//double [] sds = Utils.columnSDs(features);
		//double [] q1 = Utils.columnQ1(features);
		//double [] q3 = Utils.columnQ3(features);
		
		for (Integer id:clusteredFeatures.keySet()) {
			Cluster c = clusteredFeatures.get(id);
			System.out.println("Cluster "+id+" abundance: "+c.nElements) ;
			double [] centroid = c.calcCentroid();
			int H=0;int M =0;int L=0;
			for (int j=0;j<centroid.length;j++) {
				
				double thr0 = -1;
				
				double thr1 = thresholdPerFeature[j]; //0.3;
				System.out.println("\tF"+j+"="+centroid[j]);
				if (centroid[j]>= thr1) {
					H++;
				}else if (centroid[j]<= thr0) {
					L++;
				}else {
					M++;
				}
			}
			String centroidInterpretation = " ";
			if (H>L && H>M) {
				centroidInterpretation = HIGH_EP;
				highriskclusterfound=true;
			}
			else if (L>H && L>M)
				centroidInterpretation = " ";
			
			centroid_interpretations.put(id, centroidInterpretation);
			centroids.put(id, centroid);
			System.out.println("Centroid with ID "+id+": interpreted as "+centroidInterpretation);
		}
		
	}
	
	/*
	public File toLab() throws Exception{
		double times[] = SignalProcessing.featureTimesInSec(window4Analysis,audio);
		int nfeatures = vectorID2ClusterID.keySet().size();
		String labels [] = new String[times.length];
		for (int i=0;i<nfeatures;i++) {
			Integer id = i;
			Integer clusterID = vectorID2ClusterID.get(id);
			String interpretation = centroid_interpretations.get(clusterID);
			labels[i] = interpretation; 
		}
		File outputFile = new File(audio.getAbsolutePath().replace(".wav", ".lab"));
		System.out.println("Clustering annotation written to "+outputFile.getAbsolutePath());
		CorpusCleaner.vector2LabFile(labels, times,outputFile);
		return outputFile;
		
	}
	*/
	
	public File toLab(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec) throws Exception{
		
		//toLab(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec, String [] labels) throws Exception{
		times = SignalProcessing.featureTimesInSec(windowLengthSec,samplingFrequency,signalLength);
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
		
		if(outputFile!=null)
			CorpusCleaner.vector2LabFile(labels, times,outputFile);
		
		return outputFile;
		
	}
	
	public List<double[]> getTimes(int samplingFrequency, int signalLength, double windowLengthSec) throws Exception{
		
		toLab(null,samplingFrequency,  signalLength,  windowLengthSec);
		List<double[]> timeList = new ArrayList<double[]>();
		for (int i=0;i<(times.length-1);i++) {
			double t0 = times[i];
			double t1 = times[i+1];
			String label = labels[i];
			label = label.trim();
			if (label.length()>0) {
				double[] timeInt = {t0,t1};
				timeList.add(timeInt);
			}
		}
		return timeList;
	}
	
	public static void toLab(File outputFile, List<double[]> times) throws Exception{
		List<Double> timesN = new ArrayList<>();
		List<String> labelsN = new ArrayList<>();
		timesN.add(0d);
		labelsN.add("");
		  
		for (double [] tint : times) {
			  
			  double t0 = tint[0];
			  double t1 = tint[1];
			  String l = "A";
			  timesN.add(t0);
			  labelsN.add(l);
			  timesN.add(t1);
			  labelsN.add("");
			  
		  }
		  
		  double[] classificationTimeRanges = new double[timesN.size()];
		  for (Integer h = 0; h <classificationTimeRanges.length;h++) {
			  classificationTimeRanges[h] = timesN.get(h);
		  }

		  String[] classificationLabels = new String[labelsN.size()];

		  for (Integer h = 0; h <classificationLabels.length;h++) {
			  classificationLabels[h] = labelsN.get(h);
		  }
		  
		  CorpusCleaner.vector2LabFile(classificationLabels, classificationTimeRanges, outputFile);
	}

	public File toLabCTC(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec, double minLabelTimeSec) throws Exception{
		
		//toLab(File outputFile, int samplingFrequency, int signalLength, double windowLengthSec, String [] labels) throws Exception{
		times = SignalProcessing.featureTimesInSec(windowLengthSec,samplingFrequency,signalLength);
		double lastTime = times[times.length-1];
		
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
		highriskclusterfound=false;
		while(i<ntimes) {
			Double currTime = timesNew.get(i);
			Double nextTime = lastTime;
			if (i<(ntimes-1)) {
				nextTime = timesNew.get(i+1);
			}
			if ((nextTime-currTime)<minLabelTimeSec) {
				if (labelsNew.get(i).equals(HIGH_EP)) {
					System.out.println("WARNING: RANGE TOO SHORT: "+(nextTime-currTime)+"s");
				}
				labelsNew.set(i, " ");
			}else {
				if (labelsNew.get(i).equals(HIGH_EP)) {
					highriskclusterfound=true;
					System.out.println("OK: SUITABLE RANGE FOUND: "+(nextTime-currTime)+"s");
				}
			}
			
			i++;
		}
		
		times = new double[timesNew.size()];
		labels = new String[timesNew.size()];
        for (int j = 0; j < times.length; j++) {
            times[j] = timesNew.get(j).doubleValue();
            labels[j] = labelsNew.get(j);
        }
        
		
		CorpusCleaner.vector2LabFile(labels, times,outputFile);
		return outputFile;
	}

}
