package it.cnr.clustering;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

public class QuickKMeans {
	public int maxIterations = 10000;
	
    public static void main(String[] args) throws Exception{
        // Example usage
        double[][] data = {
            {1.0, 2.0},
            {5.0, 8.0},
            {1.5, 1.8},
            {8.0, 8.0},
            {1.0, 0.6},
            {9.0, 11.0}
        };

        int k = 2; // Number of clusters
        System.out.println("start");
        QuickKMeans cl = new QuickKMeans();
        /*
        double[][] centroids = cl.kMeans(data, k, new File("./"));

        System.out.println("Final centroids:");
        for (double[] centroid : centroids) {
            System.out.println(Arrays.toString(centroid));
        }
        */
    }

    public File kMeans(double[][] data, int k, File outputFolder) throws Exception{
        int numVectors = data.length;
        int vectorSize = data[0].length;
        System.out.println("Clustering a matrix with dimension "+numVectors+" X "+vectorSize+" over "+k+" clusters");
        // Initialize centroids randomly
        double[][] centroids = initializeCentroids(data, k);

        // Assign vectors to clusters and update centroids until convergence
        int [] assignments = new int[numVectors];
        boolean converged = false;
        int iterations = 0;
        while (!converged && iterations < maxIterations) {
            // Assign vectors to clusters
            assignments = assignToClusters(data, centroids);

            // Update centroids
            double[][] newCentroids = updateCentroids(data, assignments, k, vectorSize);

            // Check for convergence
            converged = centroids.equals(newCentroids);

            // Update centroids for the next iteration
            centroids = newCentroids;
            iterations++;
        }
        
        System.out.println("clustering converged? "+converged);
        File outputFile = new File(outputFolder, "clustering.csv");
        String columnsNames = "id,label,cluster_id,is_an_outlier\n";
		
		StringBuffer bufferRows = new StringBuffer();
		bufferRows.append(columnsNames);
		
		for (int g = 0; g < numVectors; g++) {
			int cindex = assignments[g];
			boolean isoutlier = false;
			bufferRows.append((g + 1) + ",F" + (g + 1) + "," + cindex + "," + isoutlier + "\n");
		}
		
		BufferedWriter bwx = new BufferedWriter(new FileWriter(outputFile));
		bwx.write(bufferRows.toString());
		bwx.close();
        
        
        return outputFile;
    }

    private static double[][] initializeCentroids(double[][] data, int k) {
        double[][] centroids = new double[k][];
        
        Random random = new Random();

        for (int i = 0; i < k; i++) {
            int randomIndex = random.nextInt(data.length);
            //centroids[i] = (Arrays.copyOf(data[randomIndex], data[randomIndex].length));
            centroids[i] = data[i]; //(Arrays.copyOf(data[randomIndex], data[randomIndex].length));
        }

        return centroids;
    }

    private static int [] assignToClusters(double[][] data, double[][] centroids) {
        //List<Integer> assignments = new ArrayList<>();
        int [] assignments = new int[data.length];
        int i = 0;
        for (double[] vector : data) {
            int closestCentroid = findClosestCentroid(vector, centroids);
            assignments[i] = closestCentroid;
            i++;
        }

        return assignments;
    }

    private static int findClosestCentroid(double[] vector, double[][] centroids) {
        double minDistance = Double.MAX_VALUE;
        int closestCentroid = -1;

        for (int i = 0; i < centroids.length; i++) {
            double distance = calculateDistance(vector, centroids[i]);
            if (distance < minDistance) {
                minDistance = distance;
                closestCentroid = i;
            }
        }

        return closestCentroid;
    }

    private static double calculateDistance(double[] vector1, double[] vector2) {
        double sum = 0;
        for (int i = 0; i < vector1.length; i++) {
            double diff = vector1[i] - vector2[i];
            sum += diff * diff;
        }
        return Math.sqrt(sum);
    }

    private static double[][] updateCentroids(double[][] data, int[] assignments, int k, int vectorSize) {
        
    	double[] [] newCentroids = new double[k][vectorSize];
        int[] clusterSizes = new int[k];
        double[][] clusterSums = new double[k][vectorSize];

        // Calculate sums and sizes for each cluster
        for (int i = 0; i < data.length; i++) {
            int clusterIndex = assignments[i];
            clusterSizes[clusterIndex]++;
            for (int j = 0; j < vectorSize; j++) {
                clusterSums[clusterIndex][j] += data[i][j];
            }
        }

        // Update centroids
        for (int i = 0; i < k; i++) {
            if (clusterSizes[i] > 0) {
                double[] newCentroid = new double[vectorSize];
                for (int j = 0; j < vectorSize; j++) {
                    newCentroid[j] = clusterSums[i][j] / (double)clusterSizes[i];
                }
                newCentroids[i] = newCentroid;
            }
        }

        return newCentroids;
    }
}

