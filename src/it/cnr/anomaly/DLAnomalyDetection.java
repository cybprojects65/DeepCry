package it.cnr.anomaly;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import it.cnr.clustering.Cluster;
import it.cnr.clustering.QuickKMeans;
import it.cnr.models.lstm.AnomalyDetectionVariationalAutoEncoder;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.AnomalousCryConfiguration;
import it.cnr.workflow.coroetal2024.staging.AnomalousCryDetector;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class DLAnomalyDetection extends MSEnergyEntropyAnomalyDetection {

	public DLAnomalyDetection(File audio) {
		super(audio);
	}

	boolean[] isGoodFeature = null;
	double[] scores = null;
	int[] clusters = null;
	int [] anomalousClusters = {0};
	AnomalyDetectionVariationalAutoEncoder dae;
	
	private double getScore(double[] ds, int idx) {
		return dae.final_scores[idx];
	}
	
	private int getCluster(double score, int idx) {
		//if (score<0.0038 && score>0.0037)
		int cl = clustersOfVAE[idx];
		if (cl==orderedClusters[0])
			return 0;
		else
			return -1;
	}
	
	int clustersOfVAE[];
	int orderedClusters[];
	
	private void clusterScores(double [] scores) throws Exception{
		
		double fm [][]= new double [scores.length][1];
		for (int i=0;i<fm.length;i++) {
			fm[i][0]=scores[i];
		}
		
		QuickKMeans kmeans = new QuickKMeans();
		File kmfile = new File("clustering.csv");
		kmeans.kMeans(fm, AnomalousCryConfiguration.nClusters4AnomalyDetection, kmfile.getParentFile());
		
		clustersOfVAE = new int[scores.length];
				
		List<String> km = Files.readAllLines(kmfile.toPath());
		int iii = 0;
		HashMap<Integer,Double> centroids = new HashMap<Integer, Double>();
		HashMap<Integer,Integer> centroidSizes = new HashMap<Integer, Integer>();
		
		for (String kk:km) {
			if(iii>0) {
				String [] els = kk.split(",");
				int id = Integer.parseInt(els[0])-1;
				int clid = Integer.parseInt(els[2]);
				clustersOfVAE [id] = clid;
				double score = scores[id];
				
				Double centroid = centroids.get(clid);
				Integer centroidSize = 0;
				if (centroid==null) {
					centroid = score;
					centroidSize = 1;
				}else {
					centroid = centroid+score;
					centroidSize = centroidSizes.get(clid)+1;
				}
				centroids.put(clid, centroid);
				centroidSizes.put(clid, centroidSize);
			}
			iii++;
		}
		
		List<Integer> orderedCentroidIDs = new ArrayList<>();
		List<Double> orderedCentroids = new ArrayList<>();
		List<Integer> orderedCentroidSizes = new ArrayList<>();
		for (int kkk=0;kkk<centroids.size();kkk++) {
			double avg = centroids.get(kkk);
			int size = centroidSizes.get(kkk);
			boolean found = false;
			for (int i=0;i<orderedCentroids.size();i++) {
				if (orderedCentroids.get(i)<avg) {
					orderedCentroids.add(i,avg);
					orderedCentroidIDs.add(i,kkk);
					orderedCentroidSizes.add(i,size);
					found = true;
					break;
				}
			}
			if(!found) {
				orderedCentroids.add(avg);
				orderedCentroidIDs.add(kkk);
				orderedCentroidSizes.add(size);
			}
		}
		
		orderedClusters = new int[orderedCentroidIDs.size()];
		
		for (int kkk=0;kkk<orderedCentroidIDs.size();kkk++) {
			double avg = orderedCentroids.get(kkk);
			int size = orderedCentroidSizes.get(kkk);
			int id = orderedCentroidIDs.get(kkk);
			orderedClusters[kkk] = id;
			System.out.println("C:"+id+" AVG:"+avg+" SIZE:"+size);
		}

	}
	
	private double[][] modelData(double[][] features) {
		double[] q3 = UtilsVectorMatrix.columnQ3(features);
		double[][] binaryMatrix = UtilsVectorMatrix.binariseMatrix(features, q3);
		double[][] goodMatrix = selectGoodMS(binaryMatrix);
		//train BZM
		dae = new AnomalyDetectionVariationalAutoEncoder();//new AnomalyDetectionAutoEncoder();
		try {
			System.out.println("Training ...");
			dae.train(goodMatrix, AnomalousCryConfiguration.nhidden, AnomalousCryConfiguration.nEpochs);
			System.out.println("Testing anomalies ...");
			dae.test(goodMatrix,AnomalousCryConfiguration.reconstructionNumSamples);
			System.out.println("Clustering anomalies ...");
			clusterScores(dae.final_scores);
			System.out.println("Clustering the anomalies done");
		} catch (Exception e) {
			e.printStackTrace();
			System.exit(0);
		}
		
		//test BZM
		scores = new double[features.length];
		clusters = new int[features.length];
		//anomalousClusters= {0};
		int bzi = 0;
		for (int i=0;i<features.length;i++) {
			if (isGoodFeature[i]) {
				scores[i] = getScore(goodMatrix[bzi],bzi);
				clusters[i] = getCluster(scores[i],bzi);
				bzi++;
			}else {
				scores[i] = -1;
				clusters[i] = -1;
			}
			
		}
		
		
		return binaryMatrix;
	}



	public boolean isgoodMS(double[] b) {

		for (double d : b) {
			if (d > 0)
				return true;
		}

		return false;

	}

	public double[][] selectGoodMS(double[][] binaryMatrix) {
		isGoodFeature = new boolean[binaryMatrix.length];

		List<double[]> selected = new ArrayList<>();

		for (int i = 0; i < binaryMatrix.length; i++) {

			double[] a = binaryMatrix[i];
			if (isgoodMS(a)) {
				selected.add(a);
				isGoodFeature[i] = true;
			}else
				isGoodFeature[i] = false;
		}

		double[][] selectedM = new double[selected.size()][binaryMatrix[0].length];

		for (int j = 0; j < selectedM.length; j++) {

			selectedM[j] = selected.get(j);

		}

		return selectedM;
	}

	public void classifyFeatures(double[][] features) throws Exception {
		System.out.println("Starting MS filtering for medium entropy energy segments");

		modelData(features);

		int nrow = features.length;
		
		int hepclid = 1;
		int lepclid = 0;
		highriskclusterfound = false;
		boolean lowriskclusterfound = false;

		centroids = new HashMap<Integer, double[]>();
		centroid_interpretations = new HashMap<Integer, String>();
		clusteredFeatures = new LinkedHashMap<Integer, Cluster>();
		vectorID2ClusterID = new LinkedHashMap<Integer, Integer>();
		
		for (int i = 0; i < nrow; i++) {

			String interpretation = " ";
			double time = ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift);

			if (isAnomaly(features[i],i)) {// 100 precision
				interpretation = "Anomalous";
				//if (time<10)
					//System.out.println("ENT:" + time + ":" + centroidIndicator+" S:"+scores[i]);
				if (!highriskclusterfound) {
					highriskclusterfound = true;
					centroid_interpretations.put(hepclid, interpretation);
					Cluster c = new Cluster(hepclid, false);
					clusteredFeatures.put(hepclid, c);
				}
				Cluster c = clusteredFeatures.get(hepclid);
				c.add(features[i], i);
				centroids.put(hepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, hepclid);
				clusteredFeatures.put(hepclid, c);
			} else {

				if (!lowriskclusterfound) {
					lowriskclusterfound = true;
					centroid_interpretations.put(lepclid, interpretation);
					Cluster c = new Cluster(lepclid, false);
					clusteredFeatures.put(lepclid, c);
				}
				Cluster c = clusteredFeatures.get(lepclid);
				c.add(features[i], i);
				centroids.put(lepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, lepclid);
				clusteredFeatures.put(lepclid, c);
			}

		}

		System.out.println("Number of simulated-clusters identified: " + clusteredFeatures.keySet().size());
		if (highriskclusterfound)
			System.out.println("At least one high-enegy cluster found! It contains "
					+ clusteredFeatures.get(hepclid).nElements + " elements");

	}

	private boolean isAnomaly(double[] ds, int idx) {
		if (!isGoodFeature[idx])
			return false;
		else {
			int cluster = clusters[idx];
			for (int c : anomalousClusters) {
				if (c==cluster)
					return true;
			}
			return false;
		}
	}

}
