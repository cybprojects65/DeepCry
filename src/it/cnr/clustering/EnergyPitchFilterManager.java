package it.cnr.clustering;

import java.io.File;
import java.util.HashMap;
import java.util.LinkedHashMap;

public class EnergyPitchFilterManager extends DetectionManager{
	
	
	public EnergyPitchFilterManager(File audio) {
		super(audio);
	}
	
	public void detectHighValuedFeatures(double [][] features, double[] thresholdPerFeature) throws Exception{
		System.out.println("Starting energy-pitch filtering");
		int nrow = features.length;
		int ncol = features[0].length;
		int hepclid = 1;
		int lepclid = 0;
		highriskclusterfound=false;
		boolean lowriskclusterfound=false;
		
		centroids = new HashMap<Integer, double[]>();
		centroid_interpretations = new HashMap<Integer, String>();
		clusteredFeatures = new LinkedHashMap<Integer, Cluster>();
		vectorID2ClusterID = new LinkedHashMap<Integer, Integer>();
		
		for (int i=0;i<nrow;i++) {
			int H=0;int M =0;
			for (int j=0;j<ncol;j++) {
				
				if (features[i][j]>=thresholdPerFeature[j])
					H++;
				else
					M++;
			}
			//save high energy-pitch vector
			if (H>M) {
				if (!highriskclusterfound) {
					highriskclusterfound = true;
					centroid_interpretations.put(hepclid, HIGH_EP);
					Cluster c = new Cluster(hepclid,false);
					clusteredFeatures.put(hepclid, c);
				}
				Cluster c = clusteredFeatures.get(hepclid);
				c.add(features[i], i);
				centroids.put(hepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, hepclid);
				clusteredFeatures.put(hepclid,c);
			}else {
				if (!lowriskclusterfound) {
					lowriskclusterfound = true;
					centroid_interpretations.put(lepclid, " ");
					Cluster c = new Cluster(lepclid,false);
					clusteredFeatures.put(lepclid, c);
				}
				
				Cluster c = clusteredFeatures.get(lepclid);
				c.add(features[i], i);
				centroids.put(lepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, lepclid);
				clusteredFeatures.put(lepclid,c);
			}
			
		}
		
		System.out.println("Number of simulated-clusters identified: "+clusteredFeatures.keySet().size());
		if (highriskclusterfound)
			System.out.println("At least one high-enegy cluster found! It contains "+clusteredFeatures.get(hepclid).nElements+" elements");
		
	}
	

}
