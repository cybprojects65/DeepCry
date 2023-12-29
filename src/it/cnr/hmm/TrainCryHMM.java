package it.cnr.hmm;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import it.cnr.deeplearning.DeepLearningManager;
import it.cnr.features.FeatureExtractor;
import it.cnr.workflow.configuration.Configuration;

public class TrainCryHMM {

	public static void main(String[] args) throws Exception{
		
		File hmmfile = new File(DeepLearningManager.referenceFolder,"cryhmm.hmm");
		
		DeepLearningManager dlo = new DeepLearningManager();
		Configuration config = new Configuration();
		
		HashMap<Integer, List<double[]>> references = dlo.getReferenceVectors(
				config.energyWindow4Analysis,  
				config.pitchWindow4Analysis,
				config.featurewindowsize,  
				config.featurewindowshift);
		
		int nFeaturesPerTime = dlo.nFeaturesPerTime;
		
		List<double[][]> allTimeSeries = new ArrayList<>(); 
	
		for (Integer i:references.keySet()) {
			
			List<double[]> features = references.get(i);
			FeatureExtractor fe = new FeatureExtractor();
			
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
