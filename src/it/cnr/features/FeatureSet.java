package it.cnr.features;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.List;

import it.cnr.speech.filters.ModulationSpectrogram;

public class FeatureSet implements Serializable{
	
	
	/**
	 * 
	 */
	private static final long serialVersionUID = 1L;
	
	public List<double []> timeIntervals;
	public List<double[][]> featuresPerTimeInterval;
	
	public FeatureSet(List<double []> timeIntervals, List<double[][]> featuresPerTimeInterval) {
		this.timeIntervals = timeIntervals;
		this.featuresPerTimeInterval = featuresPerTimeInterval;
	}
	
	
	
	public static FeatureSet fromModulationSpectrogram(double [][] modulationSpectrogram, double times[], String[] labels, double modSpecWindowSec) {
		
		List<double []> timeIntervals = new ArrayList<>(); 
		List<double[][]> featuresPerTimeInterval = new ArrayList<>();
		
		for (int i=0;i<times.length;i++) {
			String label = labels[i].trim();
			if (label.length()>0) {
				double time0 = times[i];
				double time1 = times[times.length-1];
				
							
				if (i< (times.length-1 ) )
					time1 = times[i+1];
				
				double [] timeRange = {time0,time1};
				timeIntervals.add(timeRange);
				
				int sampletime0 = ModulationSpectrogram.getMSIndex(time0, modSpecWindowSec);
				int sampletime1 = ModulationSpectrogram.getMSIndex(time1, modSpecWindowSec);
				int ncol = modulationSpectrogram[0].length;
				int nrow = sampletime1-sampletime0+1;
				double [][] submatrix = new double [nrow][ncol];
				
				for (int r = sampletime0; r<=sampletime1; r++) {
					for (int c = 0; c<ncol; c++) {
						submatrix[r-sampletime0][c] = modulationSpectrogram[r][c];
					}
				}
				
				featuresPerTimeInterval.add(submatrix);
			}
		}
		
		return new FeatureSet(timeIntervals, featuresPerTimeInterval);
	}
	
}
