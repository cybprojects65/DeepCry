package it.cnr.clustering;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;
import it.cnr.workflow.utils.UtilsVectorMatrix;

public class ModulationSpectrogramFilterManager extends ClassificationManager{
	
	public ModulationSpectrogramFilterManager(WorkflowConfiguration config, File audio) {
		super(config, audio);
	}
	
	
	public double [] BoltzmannAnomaly(double [][] features) {
		
		double[] q3 = UtilsVectorMatrix.columnQ3(features);
		
		
		return null;
	}
	public void classifyFeatures(double [][] features, File outputFolder, int minClusters, int maxClusters, double entropyThreshold, double lowestEntropyThreshold) throws Exception{
		System.out.println("Starting MS filtering for medium entropy energy segments");

		//Boltzmann anomaly detection
		
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
		
		for (int i=0;i<nrow;i++) {
			
			
			double entropy = UtilsMath.roundDecimal(SignalProcessing.calculateSpectralEntropy(features[i]),2);
			double energy = UtilsMath.roundDecimal(UtilsMath.mean(features[i]),2);
			double indicator =UtilsMath.roundDecimal(energy*entropy,2);
			String centroidIndicator = indicator+" ["+entropy+"/"+energy+"]";
			//System.out.println("ENT:"+ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift)+":"+centroidIndicator);
			String entropyInterpretation = " ";
			
			//if (indicator >lowestEntropyThreshold && indicator < entropyThreshold) {
			//if (energy> 17.65 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) { 
			if (energy> 17 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) {//optimal
			//if (energy> 17.5 && energy< 20.70 && entropy > 1.05 && entropy < 1.45) { //restricts the range
				System.out.println("ENT:"+ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift)+":"+centroidIndicator);
				System.out.println("\tANOMALOUS!");
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
			

		
	}
	

}
