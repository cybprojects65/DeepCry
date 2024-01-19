package it.cnr.clustering;

import java.io.BufferedWriter;
import java.io.ByteArrayInputStream;
import java.io.File;
import java.io.FileWriter;
import java.io.InputStream;

import com.rapidminer.RapidMiner;

import weka.core.Instances;
import weka.core.converters.CSVLoader;

public class XMeansClustering extends MultiKMeans{

	public XMeansClustering() {
		super();
		
	}
	public static boolean RapidMinerInitialised = false;
	
	public File clusterFeatures(double[][] features, File outputFolder, int minNofPointsToDefineACluster, int minClusters, int maxClusters) throws Exception{
		if (!RapidMinerInitialised) {
			System.setProperty("rapidminer.init.operators", "./config/operators.xml");
			//RapidMiner.init();
			System.out.println("Rapid Miner initialized");
			
			}
			System.out.println("Running XMEANS");
			
			File clusteringFile = new File(outputFolder,"clustering.csv");

			CSVLoader loader = new CSVLoader();
			StringBuffer sb = new StringBuffer();
			System.out.println("Clustering "+features.length+" features");
			for (int i = -1; i < features.length; i++) {
				for (int j = 0; j < features[0].length; j++) {
					if (i == -1)
						sb.append("F" + j);
					else
						sb.append(features[i][j]);
					if (j < (features[0].length - 1)) {
						sb.append(",");
					} else
						sb.append("\n");
				}
			}
			
			InputStream tis = new ByteArrayInputStream(sb.toString().getBytes("UTF-8"));
			loader.setSource(tis);
			//Note: this requires Java 8
			Instances id = loader.getDataSet();
			long ti = System.currentTimeMillis();
			
			System.out.println("XMeans: Clustering ...");
			XMeans xmeans = new XMeans();
			xmeans.setMaxIterations(100);
			xmeans.setMinNumClusters(minClusters);
			xmeans.setMaxNumClusters(maxClusters);
			xmeans.buildClusterer(id);
			System.out.println("XMEANS: ...ELAPSED CLUSTERING TIME: " + (System.currentTimeMillis() - ti));

			// do clustering
			Instances is = xmeans.getClusterCenters();
			int nClusters = is.numInstances();
			System.out.println("Estimated " + nClusters + " best clusters");
			int[] clusteringAssignments = xmeans.m_ClusterAssignments;
			
			String columnsNames = "id,label,cluster_id,is_an_outlier\n";
			
			StringBuffer bufferRows = new StringBuffer();
			bufferRows.append(columnsNames);
			int nrows = features.length;
			for (int k = 0; k < nrows; k++) {
				int cindex = clusteringAssignments[k];
				boolean isoutlier = false;
				bufferRows.append((k + 1) + ",F" + (k + 1) + "," + cindex + "," + isoutlier + "\n");
			}
			
			BufferedWriter bwx = new BufferedWriter(new FileWriter(clusteringFile));
			bwx.write(bufferRows.toString());
			bwx.close();
			return clusteringFile;

	}
	
	
}
