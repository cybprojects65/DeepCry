package it.cnr.workflow.configuration;

public class WorkflowConfiguration {

	public float maxSilence = 0.500f; // s
	public float minimumAudioLength = 5; // s

	public float energyWindow4Analysis = 0.1f;
	public float pitchWindow4Analysis = 0.1f;
	public float featurewindowsize = 0.3f;
	public float featurewindowshift = 0.1f;
	
	public int minNFeaturesInCluster = 5;
	public int nClasses = 2;
	public int nhidden = 3;
	public int minibatch = 150;
	public int nEpochs = 2;
	public int nHMMStates = 2;
	
	public boolean standardiseFeatures = true;
	
	public WorkflowConfiguration() {
		
	}
	
	public WorkflowConfiguration(float maxSilence, float minimumAudioLength, float energyWindow4Analysis,
			float pitchWindow4Analysis, float featurewindowsize, float featurewindowshift, int minNFeaturesInCluster,
			int nClasses, int nhidden, int minibatch, int nEpochs, boolean standardiseFeatures) {

		
		this.maxSilence= maxSilence;
		this.minimumAudioLength =minimumAudioLength;
		this.energyWindow4Analysis=energyWindow4Analysis;
		this.pitchWindow4Analysis=pitchWindow4Analysis;
		this.featurewindowsize=featurewindowsize;
		this.featurewindowshift=featurewindowshift;
		this.minNFeaturesInCluster=minNFeaturesInCluster;
		this.nClasses=nClasses;
		this.nhidden=nhidden;
		this.minibatch=minibatch;
		this.nEpochs=nEpochs;
		this.standardiseFeatures=standardiseFeatures;
		
	}
	
	
	public String toString() {
		String report =
				"#####SEGMENTATION#####\n"+
				"maxSilence="+maxSilence+"s\n"+
				"minimumAudioLength="+minimumAudioLength+"s\n"+
				"#####FEATURE EXTRACTION#####\n"+
				"energyWindow4Analysis="+energyWindow4Analysis+"s\n"+
				"pitchWindow4Analysis="+pitchWindow4Analysis+"s\n"+
				"featurewindowsize="+featurewindowsize+"s\n"+
				"featurewindowshift="+featurewindowshift+"s\n"+
				"#####CLUSTERING#####\n"+
				"minNFeaturesInCluster="+minNFeaturesInCluster+"\n"+
				"#####DEEP LEARNING#####\n"+
				"nClasses="+nClasses+"\n"+
				"nhidden="+nhidden+"\n"+
				"minibatch="+minibatch+"\n"+
				"nEpochs="+nEpochs+"\n"+
				"standardiseFeatures="+standardiseFeatures;
		return report;
		
	}
	
	public WorkflowConfiguration clone() {
		
		return new WorkflowConfiguration( maxSilence, minimumAudioLength, energyWindow4Analysis,
				pitchWindow4Analysis, featurewindowsize, featurewindowshift, minNFeaturesInCluster,
				nClasses, nhidden, minibatch, nEpochs,standardiseFeatures);
		
	}

}
