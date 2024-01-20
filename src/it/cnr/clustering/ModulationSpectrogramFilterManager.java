package it.cnr.clustering;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;
import java.util.List;

import it.cnr.models.lstm.DichotomicBoltzmannMachine;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.staging.AnomalousCryDetector;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;
import it.cnr.workflow.utils.UtilsObjects;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class ModulationSpectrogramFilterManager extends ClassificationManager{
	
	public ModulationSpectrogramFilterManager(WorkflowConfiguration config, File audio) {
		super(config, audio);
	}
	
	int nEpochs = 1000;
	int encoderNeurons = 10; //fixed
	int decoderNeurons = 5; //fixed
	int representationalNeurons = 8; //fixed
	int reconstructionNumSamples = 16; //fixed
	double minPercAnomaly = 10d;
	
	public DichotomicBoltzmannMachine boltzmannAnomaly(double [][] features) throws Exception{
		
		DichotomicBoltzmannMachine boltzmann = null;
		double[] q3 = null;
		File msquantilesFile = new File("msquantiles.bin");
		File boltzmannFile = new File("boltzmann.bin");
		
		if ( (!msquantilesFile.exists() && !boltzmannFile.exists()) || !AnomalousCryDetector.skippreprocessing) {
			q3 = UtilsVectorMatrix.columnQ3(features);
			UtilsObjects.saveObject(msquantilesFile, q3);
			
			double [][] binaryMatrix = UtilsVectorMatrix.binariseMatrix(features, q3);
			boolean isBinaryData = true;
			//TODO: CACHE THE BZM AND THE QUANTILES FOR BINARISATION!
			boltzmann = new DichotomicBoltzmannMachine();
			boltzmann.detectAnomalyComplex(binaryMatrix,minPercAnomaly, isBinaryData, encoderNeurons, decoderNeurons, representationalNeurons, nEpochs, reconstructionNumSamples);
			//System.out.println("Anomalies: " + Arrays.toString(anomalies));
			UtilsObjects.saveObject(boltzmannFile, boltzmann);
		}else {
			q3 = (double  [])UtilsObjects.loadObject(msquantilesFile);
			boltzmann = (DichotomicBoltzmannMachine) UtilsObjects.loadObject(boltzmannFile);
		}
		return boltzmann;
	}
	
	public void classifyFeatures(double [][] features, File outputFolder, int minClusters, int maxClusters, double entropyThreshold, double lowestEntropyThreshold) throws Exception{
		System.out.println("Starting MS filtering for medium entropy energy segments");

		//Boltzmann anomaly detection
		System.out.println("START OF BOLTZMANN MODEL");
		DichotomicBoltzmannMachine boltzmann = boltzmannAnomaly(features);
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
		double minEnergy = averageEnergy-0.8*averageEnergy;
		double maxEnergy = averageEnergy+0.1*averageEnergy;
		
		System.out.println("Average Energy "+averageEnergy);
		for (int i=0;i<nrow;i++) {

			double entropy = UtilsMath.roundDecimal(SignalProcessing.calculateSpectralEntropy(features[i]),2);
			double energy = UtilsMath.roundDecimal(UtilsMath.mean(features[i]),2);
			double indicator =UtilsMath.roundDecimal(energy*entropy,2);
			String centroidIndicator = indicator+" ["+entropy+"/"+energy+"]";
			//System.out.println("ENT:"+ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift)+":"+centroidIndicator);
			String entropyInterpretation = " ";
			
			int isnormal = boltzmann.isOutlier[i];
			
			//if (indicator >lowestEntropyThreshold && indicator < entropyThreshold) {//other tests
			//if (energy> 17.65 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) { //other tests
			//if (energy> 17.5 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) { //restricts the range
			//TODO: test boltzmann score around the mean anomaly: 9.2:
			/*Mean anomaly :10.62 [3.005871649465976, 8.905463806856037, 17.272090319636824]
					Mean normal :55.39 [10.259068240788212, 85.51765513680309, 93.05224798364934]*/
			if (energy> 17 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) {//optimal
			//if ( boltzmann.final_scores[i]>0 && boltzmann.final_scores[i]<10) {
			//if ( boltzmann.final_scores[i]>2) {
				//	&& energy>18
					
			//if (boltzmann.final_scores[i]>1.6 && boltzmann.final_scores[i]<10.8 
					//&& energy> 17					) {
				
				double time = ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift);
				
				System.out.println("ENT:"+time+":"+centroidIndicator);
				System.out.println("Boltzmann score: "+boltzmann.final_scores[i]+" ->A="+boltzmann.isOutlier[i]);
					
				System.out.println("\tANOMALOUS!");
				
				if (time<24)
					bolzmannScoreAnomalies.add(boltzmann.final_scores[i]);
				else
					bolzmannScoreNormal.add(boltzmann.final_scores[i]);
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
				bolzmannScoreNormal.add(boltzmann.final_scores[i]);
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
	

}
