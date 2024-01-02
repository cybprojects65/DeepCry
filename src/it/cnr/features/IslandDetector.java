package it.cnr.features;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import it.cnr.clustering.ClusteringManager;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.workflow.configuration.Configuration;

public class IslandDetector {

	ClusteringManager clustering;
	Configuration config;
	File audioFile ; 
	public IslandDetector(ClusteringManager clustering, Configuration config, File audioFile) {
		
		this.clustering = clustering;
		this.config = config;
		this.audioFile = audioFile;
	}
	
	//detects islands of high energy around a cluster of high energy-pitch 
	public File detectIslands() throws Exception{
		AudioBits bits = new AudioBits(audioFile);
		short[] signal = bits.getShortVectorAudio();
		float sfrequency = bits.getAudioFormat().getSampleRate();
		bits.ais.close();
		
		FeatureExtractor fe = new FeatureExtractor(config);
		double [] energyCurve = fe.getEnergyFeatures(audioFile);
		int nfeatures = clustering.vectorID2ClusterID.keySet().size();
		
		System.out.println("Retrieving high pich-energy cells");
		List<int[]> islands = new ArrayList<>();
		for (int i=0;i<nfeatures;i++) {
			Integer id = i;
			Integer clusterID = clustering.vectorID2ClusterID.get(id);
			String interpretation = clustering.centroid_interpretations.get(clusterID);
			if (interpretation.equals(ClusteringManager.HIGH_EP)) {
				//extend backwards
				double currentEnergy = energyCurve[i];
				int stopback = i;
				for (int k=i-1;k>=0;k--) {
					double compare_en= energyCurve[k];
					double diff = (currentEnergy-compare_en)/(compare_en);
					//if (diff < 0 || diff<0.5) {
					if (diff > 0 && diff<0.5) {
						//possibly revise: energy increase should be kept!
						stopback = k;
						break;
					}  
				}
				//extend forward
				int stopforw = i;
				for (int k=i+1;k<energyCurve.length;k++) {
					double compare_en= energyCurve[k];
					double diff = (currentEnergy-compare_en)/(compare_en);
					if (diff > 0 && diff<0.5) {
						stopforw = k;
						break;
					}  
				}
				//if we found only one vector extend the range a bit
				if (stopback == i && stopforw ==i) {
					if (i>0)
						stopback = i-1;
					if (i<(nfeatures-1))
						stopforw = i+1;
				}
				//System.out.println("Island of cell n. "+id+" has been extended of ["+(i-stopback)+";"+(stopforw-i)+"]");
				int island [] = {stopback, stopforw};
				islands.add(island);
			} 
		}
		System.out.println("N. of Captured islands "+islands.size());
		
		List<int[]> non_over_islands = new ArrayList<>();
		for (int[] island: islands) {
			
			int a = island[0];
			int b = island[1];
			
			int islA = fallsIn(a, non_over_islands);
			int islB = fallsIn(b, non_over_islands);
			if (islA==-1 && islB==-1) {
				int islandAB [] = {a,b};
				int inclusion = includes(islandAB,non_over_islands);
				if (inclusion>-1) {
					islA = inclusion;
					islB = inclusion;
				}
			}
			
			
			if (islA>-1 && islB>-1) {
				int a0 = non_over_islands.get(islA)[0];
				int b0 = non_over_islands.get(islA)[1];
				int a1 = non_over_islands.get(islB)[0];
				int b1 = non_over_islands.get(islB)[1];
				
				if (islA!=islB) {
					
					int ax = Math.min(a0,a1);
					int bx = Math.max(b0,b1);
					int islandX [] = {ax,bx};
					non_over_islands.set(islA, islandX);
					non_over_islands.remove(islB);
					non_over_islands=cleanUp(islandX, non_over_islands);
				}else {
					int ax = Math.min(a0,a);
					int bx = Math.max(b,b0);
					int islandX [] = {ax,bx};
					non_over_islands.set(islA, islandX);
					non_over_islands=cleanUp(islandX, non_over_islands);
				}
			}else if (islA>-1){
				int a0 = non_over_islands.get(islA)[0];
				int b0 = non_over_islands.get(islA)[1];
				int ax = Math.min(a0,a);
				int bx = Math.max(b,b0);
				int islandX [] = {ax,bx};
				non_over_islands.set(islA, islandX);
				non_over_islands=cleanUp(islandX, non_over_islands);
			}else if (islB>-1){
				int a0 = non_over_islands.get(islB)[0];
				int b0 = non_over_islands.get(islB)[1];
				int ax = Math.min(a0,a);
				int bx = Math.max(b,b0);
				int islandX [] = {ax,bx};
				non_over_islands.set(islB, islandX);
				non_over_islands=cleanUp(islandX, non_over_islands);
			}else {
				non_over_islands.add(island);
			}
			
		}
		
		
		islands = non_over_islands;
		System.out.println("Unifying the signals");
		List<short[]> signalIslands = new ArrayList<>();
		int totalSamples = 0;
		int idxIsl = 0;
		for (int[] island:islands) {
			int isl0 = island[0];
			int isl1 = island[1];
			
			if (isl0<isl1) {
				double t0 = (isl0*(double)config.energyWindow4Analysis);
				double t1 = (isl1*(double)config.energyWindow4Analysis);
				int i0 = Utils.timeToSamples(t0, sfrequency);
				int i1 = Utils.timeToSamples(t1, sfrequency);
				System.out.println("Signal segment n. "+idxIsl+" from "+t0+"s to "+t1+"s");
				
				short[] subsignal = Arrays.copyOfRange(signal, i0, i1+1);
				totalSamples=totalSamples+subsignal.length;
				signalIslands.add(subsignal);
				idxIsl++;
			}
			
			
		}
		
		double reduction_perc = 100d*(double)(signal.length-totalSamples)/(double)signal.length;
		System.out.println("The original signal will be reduced of "+reduction_perc+"%");
		short[] unifiedSignal = new short[totalSamples];

		int idx = 0;
		for (int i = 0;i<signalIslands.size();i++) {
				
			short[] singlesignal = signalIslands.get(i);
			for (int j=0;j<singlesignal.length;j++) {
				unifiedSignal[idx] = singlesignal[j];
				idx = idx+1;
			}
		}
		
		File singleFile = new File(audioFile.getAbsolutePath().replace(".wav", "_islands.wav"));
		AudioWaveGenerator.generateWaveFromSamplesWithSameFormat(unifiedSignal, singleFile, new AudioBits(audioFile).getAudioFormat());
		System.out.println("Energy Islands Saved to "+singleFile.getAbsolutePath()+" - Done");
		return singleFile;
	}
	
	int fallsIn(int islandIdx, List<int[]> islands) {
		int counter = 0;
		for (int[] island:islands) {
			if ((islandIdx>=island[0] && islandIdx<=island[1]) ||
					(islandIdx==(island[1]+1)) ) {
				return (counter);
			}
			counter++;
		}
		return -1;
	}
	
	int includes(int [] refisland, List<int[]> islands) {
		int counter = 0;
		for (int[] island:islands) {
			if (refisland[0]<=island[0] && refisland[1]>=island[1])
				return (counter);
			counter++;
		}
		return -1;
	}
	
	List<int[]> cleanUp(int[] refinterval, List<int[]> islands) {
		List<int[]> clean = new ArrayList<int[]>();
		//clean.add(refinterval);
		for (int[] island:islands) {
			if (refinterval[0]<island[0] && island[1]<refinterval[1]) {
				
			}else {
				clean.add(island);
			}
			
		}
		return clean;
	}
	
}