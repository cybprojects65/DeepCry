package it.cnr.clustering;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class Cluster {

	public int ID;
	public boolean isOutlier = false;
	public List<double []> features;
	public List<Integer> featureIndex;
	public int nElements = 0;
	
	public Cluster(int ID, boolean isOutlier) {
		this.ID = ID;
		this.isOutlier = isOutlier;
		features = new ArrayList<>();
		featureIndex = new ArrayList<Integer>();
	}
	
	public void add(double[] feature, Integer index) {
		Integer bestIdx = 0;
		for(Integer id:featureIndex) {
			if(id>index) {
				break;
			}
			bestIdx++; 
		}
		
		features.add(bestIdx,feature);
		featureIndex.add(bestIdx,index);
		nElements = features.size(); 
	}
	
	
	public double[] calcCentroid() {

			double[] centroid = new double[features.get(0).length];
			int nvectsincluster = features.size();

			for (int i = 0; i < centroid.length; i++) {

				for (int j = 0; j < nvectsincluster; j++) {
					centroid[i] = centroid[i] + features.get(j)[i];
				}

				centroid[i] = centroid[i] / (double) nvectsincluster;

			}
			
			return centroid;

	}
	
	
}
