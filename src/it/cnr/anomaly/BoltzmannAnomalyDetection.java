package it.cnr.anomaly;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import it.cnr.clustering.ClassificationManager;
import it.cnr.clustering.Cluster;
import it.cnr.clustering.QuickKMeans;
import it.cnr.models.lstm.DichotomicBoltzmannMachine;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.staging.AnomalousCryDetector;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;
import it.cnr.workflow.utils.UtilsObjects;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class BoltzmannAnomalyDetection extends ClassificationManager{
	
	public BoltzmannAnomalyDetection(WorkflowConfiguration config, File audio) {
		super(config, audio);
	}
	
	int nEpochs = 1000;
	int encoderNeurons = 1; //fixed
	int decoderNeurons = 1; //fixed
	int representationalNeurons = 8; //fixed
	int reconstructionNumSamples = 16; //fixed
	double minPercAnomaly = 10d;
	double [][] binaryMatrix = null;
	int [] boltzmannBelongingCluster = null;
	int [] boltzmannBelongingClusterComplete = null;
	int anomalousCluster = -1; 
	public DichotomicBoltzmannMachine boltzmannAnomaly(double [][] features) throws Exception{
		
		DichotomicBoltzmannMachine boltzmann = null;
		double[] q3 = null;
		File msquantilesFile = new File("msquantiles.bin");
		File binaryMatrixFile = new File("binaryMatrix.bin");
		File boltzmannFile = new File("boltzmann.bin");
		
		if ( (!msquantilesFile.exists() && !boltzmannFile.exists() && !binaryMatrixFile.exists()) ) // || !AnomalousCryDetector.skippreprocessing)
		//if ( !AnomalousCryDetector.skippreprocessing)
		//if (true)
		{
			q3 = UtilsVectorMatrix.columnQ3(features);
			UtilsObjects.saveObject(msquantilesFile, q3);
			
			binaryMatrix = UtilsVectorMatrix.binariseMatrix(features, q3);
			
			UtilsObjects.saveObject(binaryMatrixFile, binaryMatrix);
			
			double [][] goodMatrix = selectGoodMS(binaryMatrix);
			boolean isBinaryData = true;
			//TODO: CACHE THE BZM AND THE QUANTILES FOR BINARISATION!
			boltzmann = new DichotomicBoltzmannMachine();
			//boltzmann.detectAnomalyComplex(binaryMatrix,minPercAnomaly, isBinaryData, encoderNeurons, decoderNeurons, representationalNeurons, nEpochs, reconstructionNumSamples);
			boltzmann.detectAnomalyComplex(goodMatrix,minPercAnomaly, isBinaryData, encoderNeurons, decoderNeurons, representationalNeurons, nEpochs, reconstructionNumSamples);
			//System.out.println("Anomalies: " + Arrays.toString(anomalies));
			UtilsObjects.saveObject(boltzmannFile, boltzmann);
			
			
		}else {
			q3 = (double  [])UtilsObjects.loadObject(msquantilesFile);
			boltzmann = (DichotomicBoltzmannMachine) UtilsObjects.loadObject(boltzmannFile);
			binaryMatrix = (double [][]) UtilsObjects.loadObject(binaryMatrixFile);
		}
		
		
		double fm [][]= new double [boltzmann.final_scores.length][1];
		for (int i=0;i<fm.length;i++) {
			fm[i][0]=boltzmann.final_scores[i];
		}
		/*
		MultiKMeans kmeans = new MultiKMeans();
			File kmfile = new File("clustermskmeans.csv");
			kmeans.cluster(fm, 2, 2, 10, kmfile);
		*/
		QuickKMeans kmeans = new QuickKMeans();
		File kmfile = new File("clustering.csv");
		kmeans.kMeans(fm, 8, kmfile.getParentFile());
		
		boltzmannBelongingCluster = new int[boltzmann.final_scores.length];
		List<String> km = Files.readAllLines(kmfile.toPath());
		int iii = 0;
		HashMap<Integer,Double> centroids = new HashMap<Integer, Double>();
		HashMap<Integer,Integer> centroidSizes = new HashMap<Integer, Integer>();
		for (String kk:km) {
			if(iii>0) {
				String [] els = kk.split(",");
				int id = Integer.parseInt(els[0])-1;
				int clid = Integer.parseInt(els[2]);
				boltzmannBelongingCluster [id] = clid;
				double score = boltzmann.final_scores[id];
				
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
		
		
		for (int kkk=0;kkk<centroids.size();kkk++) {
			
			double avg = centroids.get(kkk);
			int size = centroidSizes.get(kkk);
			System.out.println("C:"+kkk+" AVG:"+avg+" SIZE:"+size);
			if (avg>600 && avg<1000) {
				anomalousCluster = kkk; //TODO : select based on NBEST
				break;
			}
		}
		
		boltzmannBelongingClusterComplete = new int[binaryMatrix.length];
		int bbi = 0;
		for (int hhh=0; hhh<binaryMatrix.length;hhh++) {
			
			if (!isgoodMS(binaryMatrix[hhh])) {
				boltzmannBelongingClusterComplete[hhh]=-1;
			}else {
				boltzmannBelongingClusterComplete[hhh]=boltzmannBelongingCluster[bbi];
				bbi++;
			}
			
		}
		
		anomalousCluster=4;
		System.out.println("Anomaly cluster: "+anomalousCluster);
		return boltzmann;
	}
	
	public double [][] selectGoodMS(double [][] binaryMatrix){
		
		List<double []> selected = new ArrayList<>();
		
		for (int i=0;i<binaryMatrix.length;i++) {
			
			double[] a = binaryMatrix[i];
			if (isgoodMS(a)) {
				selected.add(a);
			}
		}
		
		double [][] selectedM = new double[selected.size()][binaryMatrix[0].length];
		
		for (int j=0;j<selectedM.length;j++) {
			
			selectedM[j] = selected.get(j);
			
		}
		
		return selectedM;
	}
	
	public boolean isgoodMS(double [] b) {
		
		for (double d:b) {
			if (d>0)
				return true; 
		}
		
		return false;
		
	}
	
	
	public void classifyFeatures(double [][] features, File outputFolder, int minClusters, int maxClusters, double entropyThreshold, double lowestEntropyThreshold) throws Exception{
		System.out.println("Starting MS filtering for medium entropy energy segments");

		//Boltzmann anomaly detection
		System.out.println("START OF BOLTZMANN MODEL");
		DichotomicBoltzmannMachine boltzmann = boltzmannAnomaly(features);
		/*
		DichotomicBoltzmannMachine boltzmann = new DichotomicBoltzmannMachine();
		boltzmann.final_scores = UtilsVectorMatrix.initializeVector(features.length, 0);
		boltzmann.isOutlier =new int [features.length];
		*/
		
		System.out.println("END OF BOLTZMANN MODEL");
		//End of Boltzmann anomaly detection
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
		double minEnergyDetected = Double.MAX_VALUE;
		double maxEnergyDetected = 0;
		double minEntropyDetected = Double.MAX_VALUE;
		double maxEntropyDetected = 0;
		double[] spectralenergy = new double[nrow];
		double[] spectralentropy = new double[nrow];
		List<Double> bolzmannScoreAnomalies = new ArrayList<>();
		List<Double> bolzmannScoreNormal = new ArrayList<>();
		
		double totalEnergy = 0;
		for (int i=0;i<nrow;i++) {
			double energy = UtilsMath.sumVector(features[i]);
			totalEnergy+=energy;	
		}
		double averageEnergy = totalEnergy/(double) nrow;
		double minEnergy = 0.0*averageEnergy;
		double maxEnergy = averageEnergy+0.1*averageEnergy;
		
		System.out.println("Average Energy "+averageEnergy);
		int bi = 0;
		for (int i=0;i<nrow;i++) {

			double entropy = UtilsMath.roundDecimal(SignalProcessing.calculateSpectralEntropy(features[i]),2);
			double energy = UtilsMath.roundDecimal(UtilsMath.mean(features[i]),2);
			double indicator =UtilsMath.roundDecimal(energy*entropy,2);
			String centroidIndicator = indicator+" ["+entropy+"/"+energy+"]";
			//System.out.println("ENT:"+ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift)+":"+centroidIndicator);
			String entropyInterpretation = " ";
			
			//int isnormal = boltzmann.isOutlier[i];
			
			double time = ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift);
			//System.out.println("ENT:"+time+":"+centroidIndicator);
			
			
			/*Mean anomaly :10.62 [3.005871649465976, 8.905463806856037, 17.272090319636824]
					Mean normal :55.39 [10.259068240788212, 85.51765513680309, 93.05224798364934]*/
			//if (energy> 17 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) {//optimal: good precision - recall
			/*
			if (energy> 17 && energy< 20.70 && entropy > 1.05 && entropy < 1.41
					&&isgoodMS(binaryMatrix[i])
					) { //100% precision
				*/
			if (isgoodMS(binaryMatrix[i]) && 
					(isGoodCluster(boltzmannBelongingClusterComplete[i]))
							) {
			//TODO: use positive energy instead of average
			//TODO: clusterise the final scores and find correspondence between the energy/entropy ranges and the anomalies
			/*
			if (energy> minEnergy 
					&& isgoodMS(binaryMatrix[i]) 
					&& boltzmann.final_scores[bi]>74 
					&& boltzmann.final_scores[bi]<97) { //100% precision
				*/
				//System.out.println("Boltzmann score: "+boltzmann.final_scores[i]+" ->A="+boltzmann.isOutlier[i]);
				
				if (isgoodMS(binaryMatrix[i])) {
					//System.out.println("ENT:"+time+":"+centroidIndicator);
					bolzmannScoreAnomalies.add(boltzmann.final_scores[bi]);
					//System.out.println("Boltzmann score: "+boltzmann.final_scores[bi]+" ->CLUS="+boltzmannBelongingCluster[bi]+ ": ANOMAL:"+(boltzmannBelongingCluster[bi]==anomalousCluster));
					//System.out.println("\tANOMALOUS!");
					
					bi++;
				}
				entropyInterpretation = "Anomalous";
				if (!highriskclusterfound) {
					highriskclusterfound = true;
					centroid_interpretations.put(hepclid, entropyInterpretation);
					Cluster c = new Cluster(hepclid,false);
					clusteredFeatures.put(hepclid, c);
				}
				Cluster c = clusteredFeatures.get(hepclid);
				c.add(features[i], i);
				centroids.put(hepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, hepclid);
				clusteredFeatures.put(hepclid,c);
			}else {
				
				if (isgoodMS(binaryMatrix[i])) {
					bolzmannScoreNormal.add(boltzmann.final_scores[bi]);
					bi++;
				}
				
				if (!lowriskclusterfound) {
					lowriskclusterfound = true;
					centroid_interpretations.put(lepclid, entropyInterpretation);
					Cluster c = new Cluster(lepclid,false);
					clusteredFeatures.put(lepclid, c);
				}
				Cluster c = clusteredFeatures.get(lepclid);
				c.add(features[i], i);
				centroids.put(lepclid, c.calcCentroid());
				vectorID2ClusterID.put(i, lepclid);
				clusteredFeatures.put(lepclid,c);
			}
			
			if (minEnergyDetected>energy) {
				minEnergyDetected=energy;
			}if (maxEnergyDetected<energy) {
				maxEnergyDetected=energy;
			}if (minEntropyDetected>entropy) {
				minEntropyDetected=entropy;
			}if (maxEntropyDetected<entropy) {
				maxEntropyDetected=entropy;
			}
			spectralenergy[i] =energy;
			spectralentropy[i] = entropy;
			
		}
		
		System.out.println("Energy detected across the samples:["+minEnergyDetected+","+maxEnergyDetected+"]");
		System.out.println("Entropy detected across the samples:["+minEntropyDetected+","+maxEntropyDetected+"]");
		double [] qene = UtilsMath.quantiles(spectralenergy);
		double [] qent = UtilsMath.quantiles(spectralentropy);
		System.out.println("Quartiles of Energy:"+Arrays.toString(qene));
		System.out.println("Quartiles of Entropy:"+Arrays.toString(qent));
		
		System.out.println("Number of simulated-clusters identified: "+clusteredFeatures.keySet().size());
		if (highriskclusterfound)
			System.out.println("At least one high-enegy cluster found! It contains "+clusteredFeatures.get(hepclid).nElements+" elements");
			
		double bzscoresA [] = bolzmannScoreAnomalies.stream().mapToDouble(Double::doubleValue).toArray();
		double bzscoresN [] = bolzmannScoreNormal.stream().mapToDouble(Double::doubleValue).toArray();
		double quantilesA [] = UtilsMath.quantiles(bzscoresA);
		double quantilesN [] = UtilsMath.quantiles(bzscoresN);
		
		System.out.println("Meta:\nEpochs="+nEpochs+"\nencoderNeurons="+encoderNeurons+"\ndecoderNeurons="+decoderNeurons+"\nrepresentationalNeurons="+representationalNeurons+"\nreconstructionNumSamples="+reconstructionNumSamples+
				"\nminPercAnomaly="+minPercAnomaly);
		double am = UtilsMath.mean(bzscoresA);
		double an = UtilsMath.mean(bzscoresN);
		System.out.println("Mean anomaly :"+ UtilsMath.roundDecimal(am,2) +" "+Arrays.toString(quantilesA));
		System.out.println("Mean normal :"+UtilsMath.roundDecimal(an,2)  +" "+Arrays.toString(quantilesN));
		System.out.println("Diff :"+ UtilsMath.roundDecimal((an-am),2));
		
		
		
	}

	private boolean isGoodCluster(int i) {
		
		int [] goodc = {
				0,
				2,
				3,
				5,
				1,
				4
				-2
		};
		
		for (int k:goodc) {
			
			if (k==i)
				return true;
		}
		
		return false;
	}
	

}
