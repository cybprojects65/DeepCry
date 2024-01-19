package it.cnr.clustering;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.io.FileUtils;

public class MultiKMeans implements Serializable {

	/**
	 * 
	 */
	
	public MultiKMeans() {
		//KMeans.RapidMinerInitialised=false;
	}
	
	
	private static final long serialVersionUID = 1L;
	public float difformity;
	public double BIC;
	public HashMap<Integer, List<double[]>> clusters;

	public int R;
	public int K;
	public int p;
	public int M;
	public int Kstar = 0;

	public void cluster(double[][] featureList, int minElements, int minClusters, int maxClusters, File outputFile)
			throws Exception {

		HashMap<String, String> vectorLabels = new HashMap<>();
		for (int i = 0; i < featureList.length; i++) {

			vectorLabels.put("" + (i + 1), "F" + (i + 1));

		}

		double maxBIC = -Double.MAX_VALUE;
		Kstar = 0;
		double[] BICs = new double[maxClusters - minClusters + 1];
		double[] difformities = new double[maxClusters - minClusters + 1];
		if (minClusters==maxClusters)
			Kstar =minClusters;
		else {
		for (int K = minClusters; K <= maxClusters; K++) {

			//System.out.println("#######Running Kmeans clustering...K = " + K + " <= " + maxClusters);
			KMeans kmeans = new KMeans(vectorLabels);
			if (outputFile.exists())
				FileUtils.forceDelete(outputFile);

			File outfile = kmeans.compute(K, 100, 10000, minElements, featureList, outputFile);
			//System.out.println("...Done");
			int nfeatures = featureList.length;
			float uniformity = (float) nfeatures / (float) K;

			difformity = 0;
			for (String key : kmeans.pointsPerCluster.keySet()) {
				int np = kmeans.pointsPerCluster.get(key);

				difformity = difformity + (Math.abs(np - uniformity));
			}

			difformity = (float) difformity / K;

			// difformity measurement is a way to check the fairness of the kmeans
			//System.out.println("Difformity: " + difformity);
			BIC = BIC(outfile, featureList, K);
			//System.out.println("BIC: " + BIC);
			if (maxBIC < BIC) {

				maxBIC = BIC;
				Kstar = K;

			}
			BICs[K - minClusters] = BIC;
			difformities[K - minClusters] = difformity;

			//System.out.println("#######\n");
		}
		System.out.println("Search between : " + minClusters + " and " + maxClusters);
		System.out.println("Optimal BIC: " + maxBIC);
		System.out.println("Optimal K: " + Kstar);
		System.out.println("Optimal BIC: " + maxBIC);
		System.out.println("M: " + M + " R:" + R + " p: " + p);
		}
		
		KMeans kmeans = new KMeans(vectorLabels);
		kmeans.compute(Kstar, 100, 10000, minElements, featureList, outputFile);
		System.out.println(
				"BICs\n" + Arrays.toString(BICs).replace("[", "").replace("]", "").replace(",", "\t").replace(" ", ""));
		System.out.println("DIFFs\n"
				+ Arrays.toString(difformities).replace("[", "").replace("]", "").replace(",", "\t").replace(" ", ""));
		BIC = BIC(outputFile, featureList, Kstar);

	}

	public double BIC(File clusteringFile, double[][] featureList, int theoreticalK) throws Exception {

		clusters = new HashMap();

		List<String> all = Files.readAllLines(clusteringFile.toPath());
		int i = 0;
		for (String line : all) {
			if (i > 0) {

				String[] row = line.split(",");
				int id = Integer.parseInt(row[0]) - 1;
				Integer cluster = Integer.parseInt(row[2]);
				List<double[]> vecs = clusters.get(cluster);
				if (vecs == null)
					vecs = new ArrayList<>();

				vecs.add(featureList[id]);
				clusters.put(cluster, vecs);

			}

			i++;
		}

		HashMap<Integer, double[]> centroids = calcCentroids(clusters);
		Double sigma_sqr = calcSigma(clusters, centroids);
		R = R(clusters);
		K = centroids.size();
		p = p(centroids);
		M = centroids.values().iterator().next().length;
		if (K!=theoreticalK)
			return Double.MIN_VALUE;
		
		double loglike = 0;
		for (Integer key : centroids.keySet()) {

			int Rn = Rn(key, clusters);
			loglike += loglike(Rn, R, K, M, sigma_sqr);
		}

		double softneningTerm = (double) p * Math.log(R) / 2d;
		double BIC = loglike - softneningTerm;
		return BIC;
	}

	public double loglike(int Rn, int R, int K, int M, double sigma_sqr) {

		double loglike = -((double) Rn * (Math.log(Math.PI)) / 2d)
				- ((double) Rn * (double) M * Math.log(sigma_sqr) / 2d) - ((double) (Rn - K) / 2d)
				+ ((double) Rn * Math.log(Rn)) - ((double) Rn * Math.log(R));

		return loglike;
	}

	public HashMap<Integer, double[]> calcCentroids(HashMap<Integer, List<double[]>> clusters) {

		HashMap<Integer, double[]> centroids = new HashMap();

		for (Integer key : clusters.keySet()) {
			List<double[]> vecs = clusters.get(key);
			double[] centroid = new double[vecs.get(0).length];
			int nvectsincluster = vecs.size();

			for (int i = 0; i < centroid.length; i++) {

				for (int j = 0; j < nvectsincluster; j++) {
					centroid[i] = centroid[i] + vecs.get(j)[i];
				}

				centroid[i] = centroid[i] / (double) nvectsincluster;

			}

			centroids.put(key, centroid);
		}

		return centroids;

	}

	public int R(HashMap<Integer, List<double[]>> clusters) {
		int R = 0;
		for (Integer key : clusters.keySet()) {
			List<double[]> vecs = clusters.get(key);
			int Ri = vecs.size();
			R += Ri;
		}
		return R;
	}

	public int Rn(int n, HashMap<Integer, List<double[]>> clusters) {

		List<double[]> vecs = clusters.get(n);
		int Ri = vecs.size();
		return Ri;
	}

	public Double calcSigma(HashMap<Integer, List<double[]>> clusters, HashMap<Integer, double[]> centroids) {

		int K = clusters.size();
		int R = R(clusters);
		double sigma = 0;
		for (Integer key : clusters.keySet()) {
			List<double[]> vecs = clusters.get(key);

			double[] centroid = centroids.get(key);
			int nvectsincluster = vecs.size();
			for (int i = 0; i < centroid.length; i++) {

				for (int j = 0; j < nvectsincluster; j++) {
					sigma += (centroid[i] - vecs.get(j)[i]) * (centroid[i] - vecs.get(j)[i]);
				}
			}
		} // end for on keys

		sigma = sigma / (double) (R - K);

		return sigma;
	}

	public int p(HashMap<Integer, double[]> centroids) {
		int K = centroids.size();
		int M = centroids.values().iterator().next().length;
		int p = ((K - 1) + (M * K) + 1);
		return p;
	}

	public File clusterFeatures(double[][] features, File outputFolder, int minNofPointsToDefineACluster)
			throws Exception {
			return clusterFeatures(features, outputFolder, minNofPointsToDefineACluster, -1, -1);
	}
	
	public File clusterFeatures(double[][] features, File outputFolder, int minNofPointsToDefineACluster, int minClusters, int maxClusters)
			throws Exception {

		File outputFile = new File(outputFolder, "clustering.csv");
		HashMap<String, String> vectorLabels = new HashMap<>();
		for (int i = 0; i < features.length; i++) {
			vectorLabels.put("" + (i + 1), "F" + (i + 1));
		}

		System.out.println("Running Multi Kmeans clustering on " + outputFolder.getName());

		if (minClusters<0)
			minClusters = 1;
		if (maxClusters<0)
			maxClusters = features.length / minNofPointsToDefineACluster;
		
		long t0 = System.currentTimeMillis();
		System.out.println("MultiKmeans clusters " + minClusters + " - " + maxClusters + " min elements "
				+ minNofPointsToDefineACluster);

		cluster(features, minNofPointsToDefineACluster, minClusters, maxClusters, outputFile);
		long t1 = System.currentTimeMillis();
		System.out.println("...Clustering Done in " + (t1 - t0) + "ms");
		return outputFile;

	}

	public void save(File outputFile) throws Exception {

		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outputFile));
		oos.writeObject(this);
		oos.close();
		System.out.println("MultiKMeans saved");

	}

	public static MultiKMeans load(File outputFile) throws Exception {
		ObjectInputStream ois = new ObjectInputStream(new FileInputStream(outputFile));
		Object multikmeans = ois.readObject();
		ois.close();

		MultiKMeans mkm = (MultiKMeans) multikmeans;
		System.out.println("MultiKMeans loaded");

		return mkm;
	}
}
