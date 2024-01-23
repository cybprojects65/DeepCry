package it.cnr.deeplearning;

import java.io.File;
import java.io.FileWriter;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.UUID;

import it.cnr.clustering.MultiKMeans;
import it.cnr.features.EnergyPitchFeatureExtractor;
import it.cnr.hmm.HMMManager;
import it.cnr.models.lstm.DichotomicLSTM;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class DeepLearningManager {

	public HashMap<Integer, List<double[]>> goodclusters;

	public HashMap<Integer, List<double[]>> selectInterestingClusters(HashMap<Integer, List<double[]>> clusters,
			int minElements, int K) {

		//int K = clusters.size();
		
		int[] elementsPerCluster = new int[K];
		int maxElements = 0;
		for (Integer cl : clusters.keySet()) {
			
			elementsPerCluster[cl] = clusters.get(cl).size();
			if (elementsPerCluster[cl] > maxElements)
				maxElements = elementsPerCluster[cl];
		}
		System.out.println("ClusterDistribs\n" + Arrays.toString(elementsPerCluster).replace("[", "").replace("]", "")
				.replace(",", "\t").replace(" ", ""));
		System.out.println("Maximum elements per cluster=" + maxElements);
		// SELECT ONLY THOSE CLUSTERS WITH MIN ELEMENTS AND NOT FAR FROM MAXELEMENTS
		// copy old viable clusters in new clusters
		int threshold = minElements;//Math.max((int) (0.5 * maxElements), minElements);
		System.out.println("Threshold for selecting good clusters = " + threshold);
		goodclusters = new HashMap<>();
		int i = 0;
		for (int elements : elementsPerCluster) {

			if (elements >= threshold) {
				System.out.println("selected cluster #" + i + " with " + elements + " elements");
				List<double[]> clusterElements = clusters.get(i);
				goodclusters.put(i, clusterElements);
			}
			i++;
		}

		return goodclusters;

	}

	public HashMap<Integer, double[]> referenceCentroids;
	public int nFeaturesPerTime = -1;
	
	
	public static File referenceFolder = new File("./CryReferences/EnergyPitch/");
	public HashMap<Integer, List<double[]>> getReferenceVectors(float energyWindow4Analysis, float pitchWindow4Analysis,
			float featurewindowsize, float featurewindowshift) throws Exception {

		System.out.println("###RETRIEVING REFERENCE DATA###");
		File[] references = referenceFolder.listFiles();
		EnergyPitchFeatureExtractor extractor = new EnergyPitchFeatureExtractor();

		HashMap<Integer, List<double[]>> features = new HashMap<Integer, List<double[]>>();

		int i = 0;
		for (File ref : references) {

			if (ref.isDirectory()) {
				
				double[][] referenceFeatures = extractor.chunkizeTimeSeries(ref, energyWindow4Analysis,
						pitchWindow4Analysis, featurewindowsize, featurewindowshift);
				List<double[]> clust = new ArrayList<>();
				Collections.addAll(clust, referenceFeatures);

				features.put(i, clust);
			}

			i++;
		}
		nFeaturesPerTime = extractor.nFeaturesPerTime;
		MultiKMeans clusterer = new MultiKMeans();
		referenceCentroids = clusterer.calcCentroids(features);
		System.out.println("###RETRIEVED REFERENCE DATA###");
		return features;
	}

	public int calcMaxLikelyhoodClusterHMM(HashMap<Integer, List<double[]>> clusters, float energyWindow4Analysis,
			float pitchWindow4Analysis, float featurewindowsize, float featurewindowshift, int minNFeaturesInCluster, int K) throws Exception {

		getReferenceVectors(energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift);
		File hmmfile = new File(DeepLearningManager.referenceFolder,"cryhmm.hmm");
		HMMManager hmm = new HMMManager();
		hmm.loadHMM(hmmfile);
		double lowScore = -35;
		double highScore = -20;
		//double mediumScore = (highScore+lowScore)/2;
		
		HashMap<Integer, List<double[]>> features = selectInterestingClusters(clusters, minNFeaturesInCluster, K);
		double[] scores = new double[K];
		
		for (Integer clusterIdx : features.keySet()) {
			
					List<double[]> vectors = features.get(clusterIdx);
					EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
					double average = 0;
					for (double[] vec:vectors) {
						double[][] timeSeries = fe.featureToTimeSeries(vec, nFeaturesPerTime);
						double d = hmm.calcLike(timeSeries);
						average += d;
						//if (d>lowScore)
							//scores[clusterIdx]=scores[clusterIdx]+1;
					}
					
					average = average/(double)vectors.size();
					scores[clusterIdx] = average;
				}
			
		System.out.println("Cluster scores="+Arrays.toString(scores));
		
		int optimalCluster = -1;
		double optimalScore = -Double.MAX_VALUE;
		for (int i=0;i<scores.length;i++) {
			
			double score = scores[i];
			if (score!=0 && score>optimalScore)
			{
				optimalScore = score;
				optimalCluster = i;
			}
		}

		
		System.out.println("Probable cry cluster="+optimalCluster + " ("+optimalScore+")");
		
		return optimalCluster;
	}
	
	public int calcMaxLikelyhoodClusterHMM2(HashMap<Integer, List<double[]>> clusters, float energyWindow4Analysis,
			float pitchWindow4Analysis, float featurewindowsize, float featurewindowshift, int minNFeaturesInCluster, int K) throws Exception {

		getReferenceVectors(energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift);
		File hmmfile = new File(DeepLearningManager.referenceFolder,"cryhmm.hmm");
		HMMManager hmm = new HMMManager();
		hmm.loadHMM(hmmfile);
		/* 0.4
		double lowScore = -35;
		double highScore = -20;
		double minimalScore = -22;
		*/
		/*0.3
		double lowScore = -26;
		double highScore = -16; //-20; // -15.7;
		double minimalScore = -16;
		*/
		
		double lowScore = -35;
		double highScore = -20; //-20; // -15.7;
		double minimalScore = -22;
		
		//double mediumScore = (highScore+lowScore)/2;
		
		HashMap<Integer, List<double[]>> features = selectInterestingClusters(clusters, minNFeaturesInCluster, K);
		double[] scores = new double[K];
		
		for (Integer clusterIdx : features.keySet()) {
			
					List<double[]> vectors = features.get(clusterIdx);
					EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
					double average = 0;
					int ncomparisons = 0;
					for (double[] vec:vectors) {
						double[][] timeSeries = fe.featureToTimeSeries(vec, nFeaturesPerTime);
						double d = hmm.calcLike(timeSeries);
						if (d<=highScore) {
						average += d;
						//if (d>lowScore)
							//scores[clusterIdx]=scores[clusterIdx]+1;
						ncomparisons++;
						}
					}
					if (ncomparisons==0)
						average = 0;
					else
						average = average/(double)ncomparisons;
					scores[clusterIdx] = average;
				}
			
		System.out.println("Cluster scores="+Arrays.toString(scores));
		
		int optimalCluster = -1;
		double optimalScore = -Double.MAX_VALUE;
		boolean viable = false;
		for (int i=0;i<scores.length;i++) {
			
			double score = scores[i];
			if (score!=0 && score>optimalScore)
			{
				optimalScore = score;
				optimalCluster = i;
				
				if (score>minimalScore) {
					viable = true;
				}
			}
			
		}
		
		if (!viable)
			optimalCluster=-1;
		
		System.out.println("Probable cry cluster="+optimalCluster + " ("+optimalScore+")");
		
		return optimalCluster;
	}
	
	/*
	public int calcMaxLikelyhoodClusterExtensive(HashMap<Integer, List<double[]>> clusters, float energyWindow4Analysis,
			float pitchWindow4Analysis, float featurewindowsize, float featurewindowshift, int minNFeaturesInCluster) throws Exception {

		HashMap<Integer, List<double[]>> referenceVectors = getReferenceVectors(energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift);
		HashMap<Integer, List<double[]>> features = selectInterestingClusters(clusters, minNFeaturesInCluster);
		int[] scores = new int[clusters.size()];
		
		for (Integer clusterIdxRef : referenceVectors.keySet()) {
			
			List<double[]> refVectors = referenceVectors.get(clusterIdxRef);
			for (double[] ref:refVectors) {
				
				for (Integer clusterIdx : features.keySet()) {
					
					List<double[]> vectors = features.get(clusterIdx);
					
					for (double[] vec:vectors) {
						
						double angle = Utils.angle(vec, ref);
						double d = Utils.distNorm(vec, ref);
						if (angle<1)
							//System.out.println("d: "+d);
						//if (d<0.1)
							scores[clusterIdx]=scores[clusterIdx]+1;
						
					}
					
				}
				
				
			}
			
		}
			
		System.out.println("Cluster scores="+Arrays.toString(scores));
		
		int optimalCluster = -1;
		int optimalScore = 0;
		for (int i=0;i<scores.length;i++) {
			
			int score = scores[i];
			if (score>optimalScore)
			{
				optimalScore = score;
				optimalCluster = i;
			}
		}

		
		System.out.println("Probable cry cluster="+optimalCluster);
		
		return optimalCluster;
		
		
		
	}
	*/
	
	public int [] annotate (double[][] featureMatrix, File modelFolder, int nhidden,	int nClasses,	int minibatch,	int nEpochs	) throws Exception{
	
		EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
		int [] classifications = new int[featureMatrix.length];
		
		DichotomicLSTM lstm = new DichotomicLSTM();
		lstm.init(modelFolder,
				nhidden, nClasses, minibatch, nEpochs);
		
		for (int i = 0;i<featureMatrix.length;i++) {
			
			double[]v = featureMatrix[i];
			double[][] timeSeries = fe.featureToTimeSeries(v, nFeaturesPerTime);
			int c = lstm.classify(timeSeries);	
			classifications[i] = c;
			
		}
		
		return classifications;
	}
	
	public File buildTrainingSet(HashMap<Integer, List<double[]>> clusters, int trainingcluster) throws Exception{
		EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
		List<double[][]> zeroTimeSeries = new ArrayList<>();
		List<double[][]> oneTimeSeries = new ArrayList<>();
		
		for (Integer clusterIdx:clusters.keySet()) {
			
			List<double[]> vectors = clusters.get(clusterIdx);
			for (double[]v:vectors) {
				
				double[][] timeSeries = fe.featureToTimeSeries(v, nFeaturesPerTime);
				if (clusterIdx==trainingcluster) {
					oneTimeSeries.add(timeSeries);
				}else {
					zeroTimeSeries.add(timeSeries);
				}
			}
		}
		
		System.out.println("Prepared #"+oneTimeSeries.size()+" time series for 1 class");
		System.out.println("Prepared #"+zeroTimeSeries.size()+" time series for 0 class");
		
		File tempFolder = new File("temp"+UUID.randomUUID());
		
		System.out.println("Producing training and test data to "+tempFolder.getName());
		
		tempFolder.mkdir();
		File trainingFolder = new File(tempFolder,"train");
		File testFolder = new File(tempFolder,"test");
		trainingFolder.mkdir();
		testFolder.mkdir();
		File trainingFeaturesFolder = new File(trainingFolder,"features");
		trainingFeaturesFolder.mkdir();
		File trainingLabelsFolder = new File(trainingFolder,"labels");
		trainingLabelsFolder.mkdir();
		
		File testFeaturesFolder = new File(testFolder,"features");
		testFeaturesFolder.mkdir();
		File testLabelsFolder = new File(testFolder,"labels");
		testLabelsFolder.mkdir();
		
		int counter = 0;
		counter = saveTS4DL(zeroTimeSeries,0,counter,trainingFeaturesFolder,testFeaturesFolder, trainingLabelsFolder,  testLabelsFolder);
		counter = saveTS4DL(oneTimeSeries,1,counter,trainingFeaturesFolder,testFeaturesFolder, trainingLabelsFolder,  testLabelsFolder);
		System.out.println("Data preparation for DL done.");
		
		return tempFolder;
	}
	
	public double trainingAccuracy = 0;
	public String evaluation = "";
	public File trainLSTM(HashMap<Integer, List<double[]>> clusters, int trainingcluster,
			int nhidden,
	int nClasses,
	int minibatch,
	int nEpochs
	) throws Exception{
		File tempFolder = buildTrainingSet(clusters, trainingcluster); 
		DichotomicLSTM lstm = new DichotomicLSTM();
		lstm.init(tempFolder, nhidden, nClasses, minibatch, nEpochs);
		lstm.train();
		trainingAccuracy = lstm.accuracy;
		evaluation = lstm.evaluation;
		return tempFolder;
	}
	
	public int classifyLSTM(File modelFolder, double[] inputvector, int trainingcluster,int nhidden,	int nClasses,	int minibatch,	int nEpochs	) throws Exception{
		DichotomicLSTM lstm = new DichotomicLSTM();
		lstm.init(modelFolder,
				nhidden, nClasses, minibatch, nEpochs);

		//vector to TS 
		EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
		
		double[][] timeSeries = fe.featureToTimeSeries(inputvector, nFeaturesPerTime);
		
		return lstm.classify(timeSeries);
	}
	
	
	
	public int saveTS4DL(List<double[][]> timeSeries, int label, int counter, File trainingFeaturesFolder, File testFeaturesFolder, 
			File trainingLabelsFolder, File testLabelsFolder) throws Exception{
		
		for (double[][] ts:timeSeries) {
			StringBuffer sb = new StringBuffer();
			for (int i=0;i<ts.length;i++) {
				double[] row = ts[i];
				for (int j=0;j<row.length;j++) {
					double r = row[j];
					sb.append(r);
					if (j<row.length-1)
						sb.append(",");
				}
				sb.append("\n");
			}
			
			File featureFileTr = new File(trainingFeaturesFolder, counter+".csv");
			File featureFileTest = new File(testFeaturesFolder, counter+".csv");
			FileWriter fwTr = new FileWriter(featureFileTr); 
			fwTr.write(sb.toString().trim());
			fwTr.close();
			FileWriter fwTest = new FileWriter(featureFileTest); 
			fwTest.write(sb.toString().trim());
			fwTest.close();
		
			File labelFileTr = new File(trainingLabelsFolder, counter+".csv");
			File labelFileTest = new File(testLabelsFolder, counter+".csv");
			FileWriter lTr = new FileWriter(labelFileTr); 
			lTr.write(""+label);
			lTr.close();
			FileWriter lTest = new FileWriter(labelFileTest); 
			lTest.write(""+label);
			lTest.close();
			
			counter++;
		}
		
		return counter;
		
	}
	/*
	public int calcMaxLikelyhoodCluster(HashMap<Integer, List<double[]>> clusters, float energyWindow4Analysis,
			float pitchWindow4Analysis, float featurewindowsize, float featurewindowshift, int minNFeaturesInCluster) throws Exception {

		getReferenceVectors(energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize, featurewindowshift);

		HashMap<Integer, List<double[]>> features = selectInterestingClusters(clusters, minNFeaturesInCluster);
		
		MultiKMeans clusterer = new MultiKMeans();
		HashMap<Integer, double[]> centroids = clusterer.calcCentroids(features);
		HashMap<Integer, Integer> clusterScores = new HashMap<Integer, Integer>();
		
		for (Integer clusterIdxRef : referenceCentroids.keySet()) {
			double[] refcentroid = referenceCentroids.get(clusterIdxRef);
			double minDistance = Double.MAX_VALUE;
			int minIndex = -1;
			int minNFeatures = -1;
			for (Integer clusterIdx : centroids.keySet()) {

				double[] centroid = centroids.get(clusterIdx);
				double d = Utils.dist(refcentroid, centroid);
				double angle = Utils.angle(centroid, refcentroid);
				
				//if (d<minDistance) {
				if (angle<minDistance) {
					//minDistance = d;
					minDistance = angle;
					minIndex = clusterIdx;
					minNFeatures = features.get(clusterIdx).size();
				}
			}
		
			Integer score = clusterScores.get(minIndex);
			if (score == null)
				score = 1;
			else
				score = score+1;
			clusterScores.put(minIndex, score);
			
		}
		
		System.out.println("Cluster scores="+clusterScores.toString());
		
		int optimalCluster = -1;
		int optimalScore = 0;
		for (Integer clusterS : clusterScores.keySet()) {
			Integer score = clusterScores.get(clusterS);
			if (score>optimalScore)
			{
				optimalScore = score;
				optimalCluster = clusterS;
			}
		}

		
		System.out.println("Probable cry cluster="+optimalCluster);
		
		return optimalCluster;
		
		
		
	}
	*/

}
