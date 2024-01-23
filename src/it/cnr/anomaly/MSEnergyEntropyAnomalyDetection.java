package it.cnr.anomaly;

import java.io.File;
import java.util.Arrays;
import java.util.HashMap;
import java.util.LinkedHashMap;

import it.cnr.clustering.ClassificationManager;
import it.cnr.clustering.Cluster;
import it.cnr.speech.filters.ModulationSpectrogram;
import it.cnr.workflow.utils.SignalProcessing;
import it.cnr.workflow.utils.UtilsMath;

public class MSEnergyEntropyAnomalyDetection extends ClassificationManager{
	
	public MSEnergyEntropyAnomalyDetection(File audio) {
		super(audio);
	}

	public void classifyFeatures(double [][] features, double lowEnergy, double highEnergy, double lowEntropy, double highEntropy) throws Exception{
		System.out.println("Starting MS filtering for medium entropy energy segments");

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
		
		double totalEnergy = 0;
		for (int i=0;i<nrow;i++) {
			double ff [] = features[i];
			double energy = 0;
			for (double f: ff) {
				if (f>0)
					energy = energy +f;
			}
			//double energy = UtilsMath.sumVector();
			totalEnergy+=energy;	
		}
		
		double averageEnergy = totalEnergy/(double) nrow;
		double minEnergy = 0.0*averageEnergy;
		double maxEnergy = averageEnergy+0.1*averageEnergy;
		
		System.out.println("Average Energy "+averageEnergy);
		int bi = 0;
		for (int i=0;i<nrow;i++) {

			double entropy = UtilsMath.roundDecimal(SignalProcessing.calculateSpectralEntropy(features[i]),2);
			//double energy = UtilsMath.roundDecimal(UtilsMath.mean(features[i]),2);
			double energy = 0;
			for (double f: features[i]) {
				if (f>0)
					energy = energy +f;
			}
			
			double indicator =UtilsMath.roundDecimal(energy*entropy,2);
			String centroidIndicator = indicator+" ["+entropy+"/"+energy+"]";

			String entropyInterpretation = " ";
			
			double time = ModulationSpectrogram.getMSTime(i, ModulationSpectrogram.windowShift);

			//if (energy> 17 && energy< 30.70 && entropy > 1.05 && entropy < 1.45) {//optimal: good precision - recall
			/*
			if (energy> 17 && energy< 20.70 && entropy > 1.05 && entropy < 1.41
					&&isgoodMS(binaryMatrix[i])
					) { //100% precision
				*/
			//if (energy> 166 && energy<167 && entropy > 1.05 && entropy < 1.4) {//100 precision
			//if (energy> 166 && energy<167 && entropy > 0 && entropy < 1.4) {//100 precision
			if (energy> lowEnergy && energy<highEnergy && entropy > lowEntropy && entropy < highEntropy) {//100 precision
				entropyInterpretation = "Anomalous";
				System.out.println("ENT:"+time+":"+centroidIndicator);
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
			spectralenergy[i] = energy;
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
