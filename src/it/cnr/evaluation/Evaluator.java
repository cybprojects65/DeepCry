package it.cnr.evaluation;

import java.io.File;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.List;

import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.workflow.utils.SignalProcessing;

public class Evaluator {

	public void evaluate(List<File> annotations, List<File> goldenAnnotations, List<File> waveFiles) throws Exception {
		int n = waveFiles.size();

		float sfrequency = 0;

		for (int i = 0; i < n; i++) {

			File annotation = annotations.get(i);
			File goldenAnnotation = goldenAnnotations.get(i);
			if (annotation.exists() && goldenAnnotation.exists()) {
				Range[] automaticRanges = annotationToRanges(annotation);
				
				Range[] humanRanges = annotationToRanges(goldenAnnotation);

				File waveFile = waveFiles.get(i);
				AudioBits bits = new AudioBits(waveFile);
				sfrequency = bits.getAudioFormat().getSampleRate();
				int nsamples = bits.getShortVectorAudio().length;
				bits.ais.close();

				double maxTime = ((double) nsamples) / sfrequency;
				double period = 1d / sfrequency;
				double t = 0;
				while (t < maxTime) {

					boolean isAutoPositive = isInRanges(t, automaticRanges);
					boolean isHumanPositive = isInRanges(t, humanRanges);

					if (isAutoPositive && isHumanPositive)
						TP++;
					else if (isAutoPositive && !isHumanPositive) {
						// System.out.println(t);
						FP++;
					} else if (!isAutoPositive && !isHumanPositive) {
						TN++;
					} else
						FN++;

					t = t + period;
				}
			}
			// System.out.println(FP);
			// break;
		}
		
		//ADD TN from discarded files
		File folder = waveFiles.get(0).getParentFile();
		File [] allFiles = folder.listFiles();
		int pTN = TN;
		for (File f:allFiles) {
			
			if (f.getName().endsWith(".wav") && f.lastModified()==0)
			{
				
				AudioBits bits = new AudioBits(f);
				int nsamples = bits.getShortVectorAudio().length;
				bits.ais.close();
				//System.out.println("Added "+nsamples+" TN samples from discarded audio "+f.getName());
				TN +=nsamples;
			}
			
		}

		System.out.println("Added "+(TN-pTN)+" TN samples from discarded audio");
		
		calcStats(FP, FN, TP, TN);
		cohensKappa = CohensKappa(TP, TN, FN, FP);
		cohensKappaLK = kappaClassificationLandisKoch(cohensKappa);
		cohensKappaF = kappaClassificationFleiss(cohensKappa);
		System.out.println("Accuracy=" + approx(accuracy * 100));
		System.out.println("Precision=" + approx(precision * 100));
		System.out.println("Recall=" + approx(recall * 100));
		System.out.println("F1=" + approx(f1));
		System.out.println("Kappa=" + approx(cohensKappa));
		System.out.println("Kappa L&K=" + cohensKappaLK);
		System.out.println("Kappa F=" + cohensKappaF);

		System.out.println("FP=" + approx(SignalProcessing.samplesToTime(FP, sfrequency)) + "s");
		System.out.println("FN=" + approx(SignalProcessing.samplesToTime(FN, sfrequency)) + "s");
		System.out.println("TP=" + approx(SignalProcessing.samplesToTime(TP, sfrequency)) + "s");
		System.out.println("TN=" + approx(SignalProcessing.samplesToTime(TN, sfrequency)) + "s");

	}

	public double accuracy;
	public double precision;
	public double recall;
	public double f1;
	public int FP = 0;
	public int FN = 0;
	public int TP = 0;
	public int TN = 0;
	public double cohensKappa = 0;
	public String cohensKappaLK = "";
	public String cohensKappaF = "";

	public double approx(double a) {
		return (double) Math.round(a * 100) / 100d;
	}

	public void calcStats(int FP, int FN, int TP, int TN) {

		accuracy = (double) (TP + TN) / (double) (TP + FN + TN + FP);
		precision = (double) TP / (double) (TP + FP);
		recall = (double) TP / (double) (TP + FN);
		f1 = 2d * (precision * recall) / (precision + recall);

	}

	public static double CohensKappa(int TP, int TN, int FN, int FP) {
		long T = TP + FP + FN + TN;

		double Pra = (double) (TP + TN) / (double) T;
		double Pre1 = (double) (TP + FP) * (double) (TP + FN) / (double) (T * T);
		double Pre2 = (double) (TN + FN) * (double) (TN + FP) / (double) (T * T);
		double Pre = Pre1 + Pre2;
		double Kappa = (Pra - Pre) / (1d - Pre);
		return Kappa;
	}

	public boolean isInRanges(double t, Range[] ranges) {

		for (Range r : ranges) {

			if (r.isInRange(t))
				return true;

		}

		return false;
	}

	public Range[] annotationToRanges(File annotation) throws Exception {

		List<String> allLines = Files.readAllLines(annotation.toPath());
		if (allLines==null || allLines.size()==0) {
			allLines=new ArrayList<>();
			allLines.add("");
		}
		
		Range[] ranges = new Range[allLines.size()];
		int i = 0;
		for (String line : allLines) {
			if (line.trim().length()==0)
				ranges[i] = new Range(0, 0);
			else {
			String[] els = line.split(" ");
			double t0 = Double.parseDouble(els[0]);
			double t1 = Double.parseDouble(els[1]);
			ranges[i] = new Range(t0, t1);
			}
			i++;
		}

		return ranges;
	}

	public static String kappaClassificationLandisKoch(double kappa) {
		if (kappa < 0)
			return "Poor";
		else if ((kappa >= 0) && (kappa <= 0.20))
			return "Slight";
		else if ((kappa >= 0.20) && (kappa <= 0.40))
			return "Fair";
		else if ((kappa > 0.40) && (kappa <= 0.60))
			return "Moderate";
		else if ((kappa > 0.60) && (kappa <= 0.80))
			return "Substantial";
		else if (kappa > 0.80)
			return "Almost Perfect";
		else
			return "Not Applicable";
	}

	public static String kappaClassificationFleiss(double kappa) {
		if (kappa < 0)
			return "Poor";
		else if ((kappa >= 0) && (kappa <= 0.40))
			return "Marginal";
		else if ((kappa > 0.4) && (kappa <= 0.75))
			return "Good";
		else if (kappa > 0.75)
			return "Excellent";
		else
			return "Not Applicable";
	}

}
