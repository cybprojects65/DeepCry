package it.cnr.features;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import it.cnr.clustering.DetectionManager;
import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.audiofeatures.AudioWaveGenerator;
import it.cnr.workflow.configuration.AnomalousCryConfiguration;
import it.cnr.workflow.configuration.WorkflowConfiguration;
import it.cnr.workflow.coroetal2024.staging.AnomalousCryDetector;
import it.cnr.workflow.utils.SignalProcessing;

public class IslandDetector {

	DetectionManager clustering;
	//WorkflowConfiguration config;
	File audioFile ;
	private double window4Analysis; 
	
	public IslandDetector(DetectionManager clustering, WorkflowConfiguration config, File audioFile) {
		
		this.clustering = clustering;
		this.window4Analysis = config.energyWindow4Analysis;
		this.audioFile = audioFile;
	}
	
	public IslandDetector(WorkflowConfiguration config, File audioFile) {
		
		this.clustering = null;
		this.window4Analysis = config.energyWindow4Analysis;
		this.audioFile = audioFile;
	}
	
	public IslandDetector(double window4Analysis, File audioFile) {
		
		this.clustering = null;
		this.window4Analysis = window4Analysis;
		this.audioFile = audioFile;
	}

	public IslandDetector(DetectionManager clustering, double window4Analysis, File audioFile) {
		
		this.clustering = clustering;
		this.window4Analysis = window4Analysis;
		this.audioFile = audioFile;
	}


	public List<double[]> extendIntervalsThroughEnergyIslands(List<double[]> times) {
		
		EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor(window4Analysis);
		double [] energyCurve = fe.getEnergyFeatures(audioFile);
		
		List<int[]> timeIndices = new ArrayList<>();
		for (double [] t:times) {
			
			int e0 = (int) Math.floor(t[0]/window4Analysis);
			int e1 = (int) Math.floor(t[1]/window4Analysis);
			int [] newInterval = extendTimeInterval(energyCurve, e0, e1);
			timeIndices.add(newInterval);
			
		}
		
		timeIndices = mergeIntervals(timeIndices);
		
		List<double []> newtimes = new ArrayList<>();
		for (int []ti:timeIndices) {
			
			double t0 = ((double) ti[0]*window4Analysis);
			double t1 = ((double) ti[1]*window4Analysis);
			double nt [] = {t0,t1};
			newtimes.add(nt);
		}
		
		return newtimes;
	}
	
	public static double minEnergyRatio = 0.2;
	public static double maxEnergyRatio = 1000;
	public int[] extendTimeInterval(double [] energyCurve, int e0, int e1) {
		//int e0 = (int) Math.floor(t0/(double)config.energyWindow4Analysis);
		//int e1 = (int) Math.floor(t0/(double)config.energyWindow4Analysis);
		
		double currentEnergy = energyCurve[e0];
		int stopback = e0;
		for (int k=e0-1;k>=0;k--) {
			double compare_en = energyCurve[k];
			//double diff = (currentEnergy-compare_en)/(currentEnergy);
			double diff = (compare_en/currentEnergy);
			//if (diff > 0 && diff<0.5) {
			if (diff < minEnergyRatio) {
				//possibly revise: energy increase should be kept!
				stopback = k;
				break;
			}
			/*
			else if (diff > maxEnergyRatio) {
				stopback = k;
				break;
			}
			*/
			
		}
		
		int stopforw = e1;
		for (int k=e1+1;k<energyCurve.length;k++) {
			double compare_en= energyCurve[k];
			//double diff = (currentEnergy-compare_en)/(compare_en);
			double diff = (compare_en/currentEnergy);
			//if (diff > 0 && diff<0.5) {
			//if (diff > 0 && diff>0.5) {
			if (diff<minEnergyRatio) {
				stopforw = k;
				if (k<(energyCurve.length-1))
					stopforw = k+1;
				break;
			//}else if (diff < 0 && diff<-2) {
			}
			/*else if (diff>maxEnergyRatio) {
				stopforw = k;
				break;
			}
			*/
			
		}
		
		
		int[] newInterval = {stopback,stopforw};

		return newInterval;
	}
	
	
	public List<int[]> mergeIntervals(List<int[]> islands) {
		
		
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
		
		return non_over_islands;
	}
	
	//detects islands of high energy around a cluster of high energy-pitch 
	public File detectIslands() throws Exception{
		AudioBits bits = new AudioBits(audioFile);
		short[] signal = bits.getShortVectorAudio();
		float sfrequency = bits.getAudioFormat().getSampleRate();
		bits.ais.close();
		
		EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor(window4Analysis);
		double [] energyCurve = fe.getEnergyFeatures(audioFile);
		int nfeatures = clustering.vectorID2ClusterID.keySet().size();
		
		System.out.println("Retrieving high pich-energy cells");
		List<int[]> islands = new ArrayList<>();
		for (int i=0;i<nfeatures;i++) {
			Integer id = i;
			Integer clusterID = clustering.vectorID2ClusterID.get(id);
			String interpretation = clustering.centroid_interpretations.get(clusterID);
			if (interpretation.equals(DetectionManager.HIGH_EP)) {
				//extend backwards
				double currentEnergy = energyCurve[i];
				int stopback = i;
				for (int k=i-1;k>=0;k--) {
					double compare_en= energyCurve[k];
					double diff = (compare_en/currentEnergy);
					//if (diff < 0 || diff<0.5) {
					if (diff < minEnergyRatio) {
						//possibly revise: energy increase should be kept!
						stopback = k;
						break;
					}
				}
				//extend forward
				int stopforw = i;
				for (int k=i+1;k<energyCurve.length;k++) {
					double compare_en= energyCurve[k];
					double diff = (compare_en/currentEnergy);
					if (diff<minEnergyRatio) {
						stopforw = k;
						if (k<(energyCurve.length-1))
							stopforw = k+1;
						break;
					//}else if (diff < 0 && diff<-2) {
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
				double t0 = ((double) isl0*window4Analysis);
				double t1 = ((double) isl1*window4Analysis);
				int i0 = SignalProcessing.timeToSamples(t0, sfrequency);
				int i1 = SignalProcessing.timeToSamples(t1, sfrequency);
				System.out.println("Signal segment n. "+idxIsl+" from "+t0+"s to "+t1+"s");
				
				short[] subsignal = Arrays.copyOfRange(signal, i0, i1+1);
				totalSamples=totalSamples+subsignal.length;
				signalIslands.add(subsignal);
				
				short silence[] = SignalProcessing.silence(AnomalousCryConfiguration.silenceSecondsToAddBetweenDetectedHighEnergyPitchSegments, new AudioBits(audioFile).getAudioFormat().getSampleRate());
				totalSamples += silence.length;
				//add silence
				signalIslands.add(silence);
				
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
