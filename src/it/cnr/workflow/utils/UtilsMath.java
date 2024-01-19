package it.cnr.workflow.utils;

import java.math.BigInteger;
import java.util.Arrays;

public class UtilsMath {

	// ###########Math functions
		public static int powerTwoApproximation(int n) {
			int y = 0;

			while (Math.pow(2, y) <= n) {
				y++;
			}
			y--; // lower approx
			return ((int) Math.pow(2, y));

		}

		public static double dist(double[] d1, double[] d2) {

			double sum = 0;

			for (int i = 0; i < d1.length; i++) {

				sum += (d1[i] - d2[i]) * (d1[i] - d2[i]);

			}

			sum = Math.sqrt(sum);

			return sum;
		}

		public static double distNorm(double[] d1, double[] d2) {

			double sum = 0;
			double dd1 = 0;
			double dd2 = 0;
			for (int i = 0; i < d1.length; i++) {

				sum += (d1[i] - d2[i]) * (d1[i] - d2[i]);
				dd1 += (d1[i]) * (d1[i]);
				dd2 += (d2[i]) * (d2[i]);
			}

			double dist = Math.sqrt(0.5 * sum / (dd1 + dd2));

			return dist;
		}

		public static double angle(double[] d1, double[] d2) {

			double dot = 0;
			double module1 = 0;
			double module2 = 0;

			for (int i = 0; i < d1.length; i++) {

				dot += (d1[i] * d2[i]);
				module1 += (d1[i] * d1[i]);
				module2 += (d2[i] * d2[i]);
			}

			module1 = Math.sqrt(module1);
			module2 = Math.sqrt(module2);

			double angle = Math.acos(dot / (module1 * module2)) * 180 / Math.PI;

			return angle;
		}

	//rounds to the xth decimal position
		public static double roundDecimal(double number, int decimalposition) {

			double n = (double) Math.round(number * Math.pow(10.00, decimalposition)) / Math.pow(10.00, decimalposition);
			return n;
		}

		// increments a percentage o mean calculation when a lot of elements are present
		public static float incrementPerc(float perc, float quantity, int N) {

			if (N == 0)
				return quantity;

			float out = 0;
			int N_plus_1 = N + 1;
			out = (float) ((perc + ((double) quantity / (double) N)) * ((double) N / (double) N_plus_1));
			return out;

		}

		// calculates mean
		public static double mean(double[] p) {
			double sum = 0; // sum of all the elements
			for (int i = 0; i < p.length; i++) {
				sum += p[i];
			}
			return sum / p.length;
		}// end method mean

		// calculates normalized derivative
		public static double[] derivative(double[] a) {
			double[] d = new double[a.length];
			double max = 1;
			if (a.length > 0) {
				for (int i = 0; i < a.length; i++) {
					double current = a[i];
					double previous = current;
					if (i > 0) {
						previous = a[i - 1];
					}
					d[i] = current - previous;
					if (Math.abs(d[i]) > max)
						max = Math.abs(d[i]);
					// System.out.println("point "+a[i]+" derivative "+d[i]);
				}

				// normalize
				for (int i = 0; i < a.length; i++) {
					d[i] = d[i] / max;
				}
			}

			return d;
		}

		// returns a list of spikes indexes
		public static boolean[] findMaxima(double[] derivative, double threshold) {
			boolean[] d = new boolean[derivative.length];

			if (d.length > 0) {
				d[0] = false;
				for (int i = 1; i < derivative.length - 1; i++) {
					if ((derivative[i] / derivative[i + 1] < 0) && derivative[i] > 0) {
//									double ratio = Math.abs((double) derivative[i]/ (double) derivative[i+1]);
//									System.out.println("RATIO "+i+" "+Math.abs(derivative[i]));
//									if ((threshold>0)&&(ratio<threshold))
						if ((threshold > 0) && (Math.abs(derivative[i]) > threshold))
							d[i] = true;
					} else
						d[i] = false;
				}
				double max = UtilsVectorMatrix.getArgMax(derivative);
				if (max == derivative[derivative.length - 1])
					d[derivative.length - 1] = true;
				else
					d[derivative.length - 1] = false;
			}

			return d;
		}

		// returns a list of spikes indexes
		public static boolean[] findSpikes(double[] derivative, double threshold) {
			boolean[] d = new boolean[derivative.length];

			if (d.length > 0) {
				d[0] = false;
				for (int i = 1; i < derivative.length - 1; i++) {
					if (derivative[i] / derivative[i + 1] < 0) {
//								double ratio = Math.abs((double) derivative[i]/ (double) derivative[i+1]);
//								System.out.println("RATIO "+i+" "+Math.abs(derivative[i]));
//								if ((threshold>0)&&(ratio<threshold))
						if ((threshold > 0) && (Math.abs(derivative[i]) > threshold))
							d[i] = true;
					} else
						d[i] = false;
				}
				d[derivative.length - 1] = false;
			}

			return d;
		}

		// returns a list of spikes indexes
		public static boolean[] findSpikes(double[] derivative) {
			return findSpikes(derivative, -1);
		}

		public static BigInteger chunk2Index(int chunkIndex, int chunkSize) {
			return BigInteger.valueOf(chunkIndex).multiply(BigInteger.valueOf(chunkSize));
		}

		public static double[] logSubdivision(double start, double end, int numberOfParts) {

			if (end <= start)
				return null;

			if (start == 0) {
				start = 0.01;
			}
			double logStart = Math.log(start);
			double logEnd = Math.log(end);
			double step = 0;
			if (numberOfParts > 0) {

				double difference = logEnd - logStart;
				step = (difference / (double) numberOfParts);

			}
//				double [] points = new double[numberOfParts+1];
			double[] linearpoints = new double[numberOfParts + 1];

			for (int i = 0; i < numberOfParts + 1; i++) {

//					points[i] = logStart+(i*step);

				linearpoints[i] = Math.exp(logStart + (i * step));
				if (linearpoints[i] < 0.011)
					linearpoints[i] = 0;
			}

			return linearpoints;
		}

		public static double cohensKappaForDichotomy(long NumOf_A1_B1, long NumOf_A1_B0, long NumOf_A0_B1,
				long NumOf_A0_B0) {
			long T = NumOf_A1_B1 + NumOf_A1_B0 + NumOf_A0_B1 + NumOf_A0_B0;

			double Pra = (double) (NumOf_A1_B1 + NumOf_A0_B0) / (double) T;
			double Pre1 = (double) (NumOf_A1_B1 + NumOf_A1_B0) * (double) (NumOf_A1_B1 + NumOf_A0_B1) / (double) (T * T);
			double Pre2 = (double) (NumOf_A0_B0 + NumOf_A0_B1) * (double) (NumOf_A0_B0 + NumOf_A1_B0) / (double) (T * T);
			double Pre = Pre1 + Pre2;
			double Kappa = (Pra - Pre) / (1d - Pre);
			return roundDecimal(Kappa, 3);
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

		public static double scalarProduct(double[] a, double[] b) {

			double sum = 0;

			for (int i = 0; i < a.length; i++) {
				if (i < b.length)
					sum = sum + a[i] * b[i];
			}

			return sum;
		}

		public static double sumVector(double[] a) {

			double sum = 0;

			for (int i = 0; i < a.length; i++) {
				sum = sum + a[i];
			}

			return sum;
		}

		public static double[] vectorialDifference(double[] a, double[] b) {

			double[] diff = new double[a.length];

			for (int i = 0; i < a.length; i++) {
				if (i < b.length)
					diff[i] = a[i] - b[i];
				else
					diff[i] = a[i];
			}

			return diff;
		}

		public static double[] vectorialAbsoluteDifference(double[] a, double[] b) {

			double[] diff = new double[a.length];

			for (int i = 0; i < a.length; i++) {
				if (i < b.length)
					diff[i] = Math.abs(a[i] - b[i]);
				else
					diff[i] = Math.abs(a[i]);
			}

			return diff;
		}
		
		// calculates the number of elements to take from a set with inverse weight
		// respect to the number of elements
		public static int calcNumOfRepresentativeElements(int numberOfElements, int minimumNumberToTake) {
			return (int) Math.max(minimumNumberToTake, numberOfElements / Math.log10(numberOfElements));
		}

		public static double[] linearInterpolation(double el1, double el2, int intervals) {

			double step = (el2 - el1) / (double) intervals;

			double[] intervalsd = new double[intervals];
			intervalsd[0] = el1;
			for (int i = 1; i < intervals - 1; i++) {
				intervalsd[i] = el1 + step * i;
			}
			intervalsd[intervals - 1] = el2;

			return intervalsd;
		}

		private static double parabol(double a, double b, double c, double x, double shift) {
			return a * (x - shift) * (x - shift) + b * (x - shift) + c;
		}

		public static double[] inverseParabol(double a, double b, double c, double y) {

			double[] ret = { (-1d * b + Math.sqrt(b * b + 4 * a * (Math.abs(y) - c))) / (2 * a),
					(-1d * b - Math.sqrt(b * b + 4 * a * (Math.abs(y) - c))) / (2 * a) };
			return ret;
		}

		public static double logaritmicTransformation(double y) {
			y = Math.abs(y);
			if (y == 0)
				return -Double.MAX_VALUE;
			else
				return Math.log10(y);
		}

		// the parabol is centered on the start Point
		public static double[] parabolicInterpolation(double startP, double endP, int intervals) {

			double start = startP;
			double end = endP;
			double shift = start;

			double a = 1000d;
			double b = 0d;
			double c = 0d;
			double parabolStart = parabol(a, b, c, start, shift);
			if (start < 0)
				parabolStart = -1 * parabolStart;

			double parabolEnd = parabol(a, b, c, end, start);
			if (end < 0)
				parabolEnd = -1 * parabolEnd;

			double step = 0;
			if (intervals > 0) {
				double difference = Math.abs(parabolEnd - parabolStart);
				step = (difference / (double) intervals);
			}

			double[] linearpoints = new double[intervals];

			linearpoints[0] = startP;
			// System.out.println("Y0: "+parabolStart);
			for (int i = 1; i < intervals - 1; i++) {
				double ypoint = 0;
				if (end > start)
					ypoint = parabolStart + (i * step);
				else
					ypoint = parabolStart - (i * step);
				// System.out.println("Y: "+ypoint);
				double res[] = inverseParabol(a, b, c, Math.abs(ypoint));
				// System.out.println("X: "+linearpoints[i]);
				if (ypoint < 0)
					linearpoints[i] = res[1] + shift;
				else
					linearpoints[i] = res[0] + shift;
			}

			linearpoints[intervals - 1] = endP;
			return linearpoints;
		}

		// distributes uniformly elements in parts
		public static int[] takeChunks(int numberOfElements, int partitionFactor) {
			int[] partitions = new int[1];
			if (partitionFactor <= 0) {
				return partitions;
			} else if (partitionFactor == 1) {
				partitions[0] = numberOfElements;
				return partitions;
			}

			int chunksize = numberOfElements / (partitionFactor);
			int rest = numberOfElements % (partitionFactor);
			if (chunksize == 0) {
				partitions = new int[numberOfElements];
				for (int i = 0; i < numberOfElements; i++) {
					partitions[i] = 1;
				}
			} else {
				partitions = new int[partitionFactor];
				for (int i = 0; i < partitionFactor; i++) {
					partitions[i] = chunksize;
				}

				for (int i = 0; i < rest; i++) {
					partitions[i]++;
				}

			}

			return partitions;
		}

		public static int chunkize(int numberOfElements, int partitionFactor) {
			int chunksize = numberOfElements / partitionFactor;
			int rest = numberOfElements % partitionFactor;
			if (chunksize == 0)
				chunksize = 1;
			else if (rest != 0)
				chunksize++;
			/*
			 * int numOfChunks = numberOfElements / chunksize; if ((numberOfElements %
			 * chunksize) != 0) numOfChunks += 1;
			 */

			return chunksize;
		}

		public static double[] uniformSampling(double min, double max, int maxElementsToTake) {
			double step = (max - min) / (double) (maxElementsToTake - 1);
			double[] samples = new double[maxElementsToTake];

			for (int i = 0; i < samples.length; i++) {
				double value = min + i * step;
				if (value > max)
					value = max;
				samples[i] = value;
			}

			return samples;
		}

		public static int[] uniformIntegerSampling(double min, double max, int maxElementsToTake) {
			double step = (max - min) / (double) (maxElementsToTake - 1);
			int[] samples = new int[maxElementsToTake];

			for (int i = 0; i < samples.length; i++) {
				double value = min + i * step;
				if (value > max)
					value = max;
				samples[i] = (int) value;
			}

			return samples;
		}

		// finds the best subdivision for a sequence of numbers
		public static double[] uniformDivide(double max, double min, double[] points) {
				int maxintervals = 10;
				int n = maxintervals;

				boolean subdivisionOK = false;
				double gap = (max - min) / n;

				// search for the best subdivision: find the best n
				while (!subdivisionOK) {
					// System.out.println("*************************");
					boolean notcontains = false;
					// take the gap interval to test
					for (int i = 0; i < n; i++) {
						double rightmost = 0;
						// for the last border take a bit more than max
						if (i == n - 1)
							rightmost = max + 0.01;
						else
							rightmost = min + gap * (i + 1);
						// if the interval doesn't contain any point discard the subdivision
						if (!UtilsVectorMatrix.intervalContainsPoints(min + gap * i, rightmost, points)) {
							notcontains = true;
							break;
						}
					}

					// if there are empty intervals and there is space for another subdivision
					// proceed
					if (notcontains && n > 0) {
						n--;
						gap = (max - min) / n;
					}
					// otherwise take the default subdivision
					else if (n == 0) {
						n = maxintervals;
						subdivisionOK = true;
					}
					// if all the intervals are non empty then exit
					else
						subdivisionOK = true;
				}

				// once the best n is found build the intervals
				double[] intervals = new double[n];
				for (int i = 0; i < n; i++) {
					if (i < n - 1)
						intervals[i] = min + gap * (i + 1);
					else
						intervals[i] = Double.POSITIVE_INFINITY;
				}

				return intervals;
			}

		public static double[] quantiles(double[] data) {

			Arrays.sort(data);
			double q1 = calculateQuartile(data, 0.25);
			double q2 = calculateQuartile(data, 0.5);
			double q3 = calculateQuartile(data, 0.75);

			double[] quants = new double[3];
			quants[0] = q1;
			quants[1] = q2;
			quants[2] = q3;
			return quants;
		}

		public static double calculateQuartile(double[] data, double percentile) {
			int n = data.length;
			double index = percentile * (n - 1) + 1;

			if (index % 1 == 0) {
				// If the index is an integer, return the corresponding element
				return data[(int) index - 1];
			} else {
				// If the index is not an integer, interpolate between two adjacent elements
				int lowerIndex = (int) Math.floor(index);
				int upperIndex = (int) Math.ceil(index);

				double lowerValue = data[lowerIndex - 1];
				double upperValue = data[upperIndex - 1];

				return lowerValue + (index - lowerIndex) * (upperValue - lowerValue);
			}
		}

		// calculates the frequency distribution for a set of points respect to a set of
		// intervals
		public static double[] calcFrequencies(double[] interval, double[] points) {
			int intervs = interval.length;
			int npoints = points.length;
			double[] frequencies = new double[intervs];
			for (int i = 0; i < intervs; i++) {

				for (int j = 0; j < npoints; j++) {

					if (((i == 0) && (points[j] < interval[i]))
							|| ((i == intervs - 1) && (points[j] >= interval[i - 1]) && (points[j] <= interval[i]))
							|| ((i > 0) && (points[j] >= interval[i - 1]) && (points[j] < interval[i]))) {
						// System.out.println("(" + (i == 0 ? "" : interval[i - 1]) + "," + interval[i]
						// + ")" + " - " + points[j]);
						frequencies[i] = frequencies[i] + 1;
					}
				}
			}

			return frequencies;
		}
}
