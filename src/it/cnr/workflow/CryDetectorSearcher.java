package it.cnr.workflow;

import java.io.File;
import java.io.FileWriter;

public class CryDetectorSearcher {

	public Configuration findOptimalConfiguration(File audio) throws Exception {

		CryDetector cryd = null;
		Configuration optimalConfig = new Configuration();
		Configuration curreConfiguration = optimalConfig.clone();

		double optimalf1 = 0;
		/*
		float maxSilences[] = { 0.500f };
		float minimumAudioLengths[] = { 5 };
		float energyWindow4Analyses[] = { 0.1f };
		float pitchWindow4Analyses[] = { 0.1f };
		float featurewindowsizes[] = { 0.3f, 0.4f, 0.5f };
		float featurewindowshifts[] = { 0.1f };
		int minNFeaturesInClusters[] = { 3 };
		int nClasses = 2;
		int nhiddens[] = { 3, 5, 10, 20};
		int minibatches[] = {5,10,20,50,100,150,200};
		int nEpochss[] = { 2, 5, 7, 10, 20, 30};
		*/
		/*
		float maxSilences[] = { 0.500f };
		float minimumAudioLengths[] = { 5, 7 };
		
		float energyWindow4Analyses[] = { 0.1f, 0.2f, 0.05f };
		float pitchWindow4Analyses[] = { 0.1f };
		
		float featurewindowsizes[] = { 0.3f };
		float featurewindowshifts[] = { 0.1f };
		int minNFeaturesInClusters[] = { 2,3,5 };
		
		int nClasses = 2;
		int nhiddens[] = { 3};
		int minibatches[] = {150};
		int nEpochss[] = { 2};
		*/
		float maxSilences[] = { 0.500f };
		float minimumAudioLengths[] = { 5};
		
		float energyWindow4Analyses[] = { 0.1f};
		float pitchWindow4Analyses[] = { 0.1f };
		
		float featurewindowsizes[] = { 0.3f };
		float featurewindowshifts[] = { 0.1f };
		int minNFeaturesInClusters[] = { 5 };
		
		int nClasses = 2;
		int nhiddens[] = { 3};
		int minibatches[] = {150};
		int nEpochss[] = { 2};
		
		long ncombinations = maxSilences.length*minimumAudioLengths.length*
				energyWindow4Analyses.length*pitchWindow4Analyses.length*
				featurewindowsizes.length*featurewindowshifts.length*minNFeaturesInClusters.length*
				nhiddens.length*minibatches.length*nEpochss.length;
		
		
		double estimatedTime = 12*ncombinations/60;
		
		System.out.println("N combinations to test="+ncombinations+" Estimated Time="+estimatedTime+" min");
		
		
		//System.exit(0);
		try {
			if (ncombinations==1) {
				curreConfiguration = new Configuration(maxSilences[0], minimumAudioLengths[0],
						energyWindow4Analyses[0], energyWindow4Analyses[0], featurewindowsizes[0],
						featurewindowshifts[0], minNFeaturesInClusters[0], nClasses, nhiddens[0],
						minibatches[0], nEpochss[0]);
				
				optimalConfig = curreConfiguration.clone();
				cryd = new CryDetector(optimalConfig);
				cryd.detect(audio);
				cryd.eval(audio);
			}
			else {	
		long idx = 0; 
		for (float maxSilence : maxSilences) {

			for (float minimumAudioLength : minimumAudioLengths) {

				for (float energyWindow4Analysis : energyWindow4Analyses) {
					float pitchWindow4Analysis = energyWindow4Analysis;

					for (float featurewindowsize : featurewindowsizes) {

						for (float featurewindowshift : featurewindowshifts) {

							for (int minNFeaturesInCluster : minNFeaturesInClusters) {
								for (int nhidden : nhiddens) {
									for (int minibatch : minibatches) {
										for (int nEpochs : nEpochss) {
											long t0 = System.currentTimeMillis();
											System.out.println("####################################################");
											curreConfiguration = new Configuration(maxSilence, minimumAudioLength,
													energyWindow4Analysis, pitchWindow4Analysis, featurewindowsize,
													featurewindowshift, minNFeaturesInCluster, nClasses, nhidden,
													minibatch, nEpochs);
											
											cryd = new CryDetector(curreConfiguration);
											cryd.detect(audio);
											cryd.eval(audio);
											if (cryd.f1 > optimalf1) {
												optimalConfig = curreConfiguration.clone();
												optimalf1 = cryd.f1;
											}
											System.out.println("Results obtained with\n"+curreConfiguration);
											System.out.println("####################################################\n\n");
											long t1 = System.currentTimeMillis();
											
											System.out.println("Elapsed: "+(t1-t0)+" ms");
													
											cryd = null;
											System.gc();
											idx++;
											
											System.out.println("STATUS="+(idx*100d/ncombinations)+" %");
											
											Thread.sleep(2000);
										} // nEpochs
									} // minibatch
								} // nhidden
							} // minNFeaturesInCluster
						} // featurewindowshift
					} // featurewindowsize
				} // energyWindow4Analysis
			} // minimumAudioLength
		} // maxsilence

			}
		
		
		System.out.println("Optimal Configuration:\n" + optimalConfig.toString());
		
		} catch(Exception e) {
			e.printStackTrace();
			System.out.println("Last configuration used\n"+curreConfiguration);
			
		}
		
		String header= "file"+","+
				"snr"+","+
				"maxSilence"+","+
				"minimumAudioLength"+","+
				"energyWindow4Analysis"+","+
				"pitchWindow4Analysis"+","+
				"featurewindowsize"+","+
				"featurewindowshift"+","+
				"minNFeaturesInCluster"+","+
				"nClasses"+","+
				"nhidden"+","+
				"minibatch"+","+
				"nEpochs"+","+
				"accuracy"+","+
				"precision"+","+
				"recall"+","+
				"f1"+","+
				"tp"+","+
				"tn"+","+
				"fp"+","+
				"fn"+","+
				"kappa"+","+
				"kappa_lk"+","+
				"kappa_f"
				+ "\n";
		
		String log =audio.getAbsolutePath()+","+
					cryd.SNR+","+
					optimalConfig.maxSilence+","+
					optimalConfig.minimumAudioLength+","+
					optimalConfig.energyWindow4Analysis+","+
					optimalConfig.pitchWindow4Analysis+","+
					optimalConfig.featurewindowsize+","+
					optimalConfig.featurewindowshift+","+
					optimalConfig.minNFeaturesInCluster+","+
					optimalConfig.nClasses+","+
					optimalConfig.nhidden+","+
					optimalConfig.minibatch+","+
					optimalConfig.nEpochs+","+
					cryd.accuracy+","+
					cryd.precision+","+
					cryd.recall+","+
					cryd.f1+","+
					cryd.evaluator.TP+","+
					cryd.evaluator.TN+","+
					cryd.evaluator.FP+","+
					cryd.evaluator.FN+","+
					cryd.evaluator.cohensKappa+","+
					cryd.evaluator.cohensKappaLK+","+
					cryd.evaluator.cohensKappaF
					;
		
		File resultSummary = new File("Res_"+audio.getName()+".csv"); 
		FileWriter fw = null;
		
		if (!resultSummary.exists()) {
			fw = new FileWriter(resultSummary);
			fw.write(header);
		}else
			fw = new FileWriter(resultSummary,true);
		
		fw.write(log+"\n");
		fw.close();
				
		return optimalConfig;
			
	}

}
