package it.cnr.hmm;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.EnergyPitchFeatureExtractor;
import it.cnr.workflow.configuration.WorkflowConfiguration;

public class TrainCryHMM {

	public static void main(String[] args) throws Exception{
		
		File hmmfile = new File(DeepLearningManager.referenceFolder,"cryhmm.hmm");
		
		DeepLearningManager dlo = new DeepLearningManager();
		WorkflowConfiguration config = new WorkflowConfiguration();
		
		HashMap<Integer, List<double[]>> references = dlo.getReferenceVectors(
				config.energyWindow4Analysis,  
				config.pitchWindow4Analysis,
				config.featurewindowsize,  
				config.featurewindowshift);
		
		int nFeaturesPerTime = dlo.nFeaturesPerTime;
		
		List<double[][]> allTimeSeries = new ArrayList<>(); 
	
		for (Integer i:references.keySet()) {
			
			List<double[]> features = references.get(i);
			EnergyPitchFeatureExtractor fe = new EnergyPitchFeatureExtractor();
			
			for (double[] f:features) {
				double[][] timeSeries = fe.featureToTimeSeries(f, nFeaturesPerTime);
				allTimeSeries.add(timeSeries);
			}
		}
		
		System.out.println("Time Series used for training="+allTimeSeries.size());
		
		HMMManager hmm = new HMMManager(config.nHMMStates);
		
		hmm.trainHMM(allTimeSeries);
		hmm.saveHMM(hmmfile);
		
		HMMManager hmm1 = new HMMManager();
		hmm1.loadHMM(hmmfile);
		
		double minScore = Double.MAX_VALUE;
		double maxScore = -Double.MAX_VALUE;
		
		for (double[][] timeseries:allTimeSeries) {
			
			double score = hmm1.calcLike(timeseries);
			System.out.println("->"+score);
			if (score>maxScore)
				maxScore = score;
			if (score<minScore)
				minScore = score;
		}
		
		
		System.out.println("min score="+minScore+" ; max score="+maxScore);
	}
}
