package it.cnr.hmm;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import be.ac.ulg.montefiore.run.jahmm.Hmm;
import be.ac.ulg.montefiore.run.jahmm.ObservationVector;
import be.ac.ulg.montefiore.run.jahmm.OpdfMultiGaussianFactory;
import be.ac.ulg.montefiore.run.jahmm.learn.KMeansLearner;

public class HMMManager {

	public int numberOfStates;
	public Hmm<ObservationVector> hmm;
	
	public List<ObservationVector> timeSeriesToObservations(double[][] timeSeries) {
		
		List<ObservationVector> sequence = new ArrayList<ObservationVector>();

		for (int i = 0; i < timeSeries.length; i++) {
			
			double vett[] = timeSeries[i];
			ObservationVector obs = new ObservationVector(vett);
			sequence.add(obs);
		}

		return sequence;
	}
	
	public List<List<ObservationVector>> timeSeriesListToObservationList(List<double[][]> timeSeriesList){
		
		List<List<ObservationVector>> allObs = new ArrayList<List<ObservationVector>>();
		
		for (double[][] timeSeries : timeSeriesList) {
			List<ObservationVector> ts = timeSeriesToObservations(timeSeries);
			allObs.add(ts);
		}
		
		return allObs;
		
	} 
	
	public HMMManager(int numberOfStates) {
		this.numberOfStates = numberOfStates;
	}
	
	public HMMManager() {
	
	}
	
	public Hmm<ObservationVector> trainHMM(List<double[][]> timeSeriesList){
		
		List<List<ObservationVector>> allObs = timeSeriesListToObservationList(timeSeriesList);
		int numFeats = allObs.get(0).get(0).dimension();
		
		
		hmm = new Hmm<ObservationVector>(numberOfStates, new OpdfMultiGaussianFactory(numFeats));
		KMeansLearner<ObservationVector> bwl = new KMeansLearner<ObservationVector>(numberOfStates, new OpdfMultiGaussianFactory(numFeats), allObs);
		hmm = bwl.learn();
		
		return hmm;
	}
	
	public double calcLike(double[][] X) {
		List<ObservationVector> oseq = timeSeriesToObservations(X);
		// apply Viterbi
		double like = hmm.lnProbability(oseq);
		return like;
	}
	
	public void saveHMM(File outputFile) throws Exception{
		ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outputFile));
		oos.writeObject(hmm);
		oos.close();
		System.out.println("HMM saved");
		
	}
	
	
	@SuppressWarnings("unchecked")
	public Hmm<ObservationVector> loadHMM(File outputFile) throws Exception{
		
		ObjectInputStream ois = new ObjectInputStream(
				new FileInputStream(outputFile));
		Object hmmobs = ois.readObject();
		
		hmm = (Hmm<ObservationVector>) hmmobs;
		numberOfStates = hmm.nbStates();
		
		ois.close();
		System.out.println("HMM loaded");
		
		return hmm;
		
	}
}
