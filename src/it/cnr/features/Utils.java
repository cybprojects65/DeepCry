package it.cnr.features;

import java.io.File;
import java.math.BigInteger;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.DataRow;
import com.rapidminer.example.table.ExampleTable;

import it.cnr.speech.audiofeatures.AudioBits;
import it.cnr.speech.clustering.BigSamplesTable;

public class Utils {

	public static double[][] subsetRows(double [][] matrix, int row0, int row1) {
		
		double [][] sub = new double [row1-row0+1][matrix.length];
		
		for (int i = row0 ; i<=row1;i++) {
			
			sub[i-row0] = matrix[i];
			
		}
		
		return sub;
		
	}
	
	public static int powerTwoApproximation(int n) {
		int y = 0;
		
		while(Math.pow(2,y)<= n) {
				y++;		
		}
		y--; //lower approx
		return ( (int) Math.pow(2, y));
		
	}
	
	public static double [] columnToVector (double [][] matrix, int i) {
		
		double [] vector = new double[matrix.length];
		for (int k=0;k<matrix.length;k++) {
			vector [k] = matrix[k][i];
		}
		return vector;
	}
	
	public static double columnMean(double [][] matrix, int i) {
		
		double [] vector = columnToVector(matrix,i);
		return mean(vector);
	}
	
	public static double columnSD(double [][] matrix, int i) {
		
		double [] vector = columnToVector(matrix,i);
		double variance = com.rapidminer.tools.math.MathFunctions.variance(vector, Double.NEGATIVE_INFINITY);
		return Math.sqrt(variance);
	}

	public static double [] columnQ1(double [][] matrix) {
		
		int ncol = matrix[0].length;
		double[]o = new double[ncol];
		
		for (int j=0;j<ncol;j++) {
			double [] vector = columnToVector(matrix,j);
			double q1 = quantiles(vector)[0];
			o[j] = q1;
		}
		return o;
	}

	public static double [] columnQ3(double [][] matrix) {
		
		int ncol = matrix[0].length;
		double[]o = new double[ncol];
		
		for (int j=0;j<ncol;j++) {
			double [] vector = columnToVector(matrix,j);
			double q3 = quantiles(vector)[2];
			o[j] = q3;
		}
		return o;
	}

	public static double [] columnMeans(double [][] matrix) {
		
		int ncol = matrix[0].length;
		double[]colmeans = new double[ncol];
		
		for (int j=0;j<ncol;j++) {
			
			double mean = columnMean(matrix, j);
			colmeans[j] = mean;
		}
		return colmeans;
	}

	public static double [] columnSDs(double [][] matrix) {
		
		int ncol = matrix[0].length;
		double[]colsds = new double[ncol];

		for (int j=0;j<ncol;j++) {
			double sd= columnSD(matrix, j);
			colsds[j] = sd;
		}
		
		return colsds;
	}

	public static double dist(double[] d1,double[] d2) {
		
		double sum = 0;
		
		for (int i=0;i<d1.length;i++) {
			
			sum += (d1[i]-d2[i])*(d1[i]-d2[i]);
			
		}
		
		sum = Math.sqrt(sum);
		
		return sum;
	}
	
	public static double samplesToTime(int samples, double fs) {
		
		return (double) samples / fs;
		
	}
	
	public static int timeToSamples(double time, double fs) {
		
		return (int) Math.round(fs * time);
		
	}

	public static double distNorm(double[] d1,double[] d2) {
		
		double sum = 0;
		double dd1 = 0;
		double dd2 = 0;
		for (int i=0;i<d1.length;i++) {
			
			sum += (d1[i]-d2[i])*(d1[i]-d2[i]);
			dd1 += (d1[i])*(d1[i]);
			dd2 += (d2[i])*(d2[i]);
		}
		
		double dist = Math.sqrt(0.5*sum/(dd1+dd2)); 		
		
		
		return dist;
	}

	public static double angle(double[] d1,double[] d2) {
		
		double dot = 0;
		double module1 = 0;
		double module2 = 0;
		
		for (int i=0;i<d1.length;i++) {
			
			dot += (d1[i]*d2[i]);
			module1 += (d1[i]*d1[i]);
			module2 += (d2[i]*d2[i]);
		}
		
		module1 = Math.sqrt(module1);
		module2 = Math.sqrt(module2);
		
		double angle = Math.acos(dot/(module1*module2))*180/Math.PI;
		
		return angle;
	}

	//rounds to the xth decimal position
		public static double roundDecimal(double number,int decimalposition){
			
			double n = (double)Math.round(number * Math.pow(10.00,decimalposition))/Math.pow(10.00,decimalposition);
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

		public static ArrayList<Integer> generateRandoms(int numberOfRandoms, int min, int max) {

			ArrayList<Integer> randomsSet = new ArrayList<Integer>();
			// if number of randoms is equal to -1 generate all numbers
			if (numberOfRandoms == -1) {
				for (int i = min; i < max; i++) {
					randomsSet.add(i);
				}
			} else {
				int numofrandstogenerate = 0;
				if (numberOfRandoms <= max) {
					numofrandstogenerate = numberOfRandoms;
				} else {
					numofrandstogenerate = max;
				}

				if (numofrandstogenerate == 0) {
					randomsSet.add(0);
				} else {
					for (int i = 0; i < numofrandstogenerate; i++) {

						int RNum = -1;
						RNum = (int) ((max) * Math.random()) + min;

						// generate random number
						while (randomsSet.contains(RNum)) {
							RNum = (int) ((max) * Math.random()) + min;
							// AnalysisLogger.getLogger().debug("generated " + RNum);
						}

						// AnalysisLogger.getLogger().debug("generated " + RNum);

						if (RNum >= 0)
							randomsSet.add(RNum);
					}

				}
			}

			return randomsSet;
		}

		public static int[] generateSequence(int elements) {
			int[] sequence = new int[elements];
			for (int i = 0; i < elements; i++) {
				sequence[i] = i;
			}
			return sequence;
		}

		public static BigInteger chunk2Index(int chunkIndex, int chunkSize) {

			return BigInteger.valueOf(chunkIndex).multiply(BigInteger.valueOf(chunkSize));

		}

		// calculates mean
		public static double mean(double[] p) {
			double sum = 0; // sum of all the elements
			for (int i = 0; i < p.length; i++) {
				sum += p[i];
			}
			return sum / p.length;
		}// end method mean

		//calculates normalized derivative
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
					if (Math.abs(d[i])>max)
						max = Math.abs(d[i]); 
					// System.out.println("point "+a[i]+" derivative "+d[i]);
				}
				
				//normalize
				for (int i = 0; i < a.length; i++) {
					d[i] = d[i]/max;
				}
			}

			return d;
		}

		// returns a list of spikes indexes
		public static boolean[] findMaxima(double[] derivative,double threshold) {
				boolean[] d = new boolean[derivative.length];

				if (d.length > 0) {
					d[0] = false;
					for (int i = 1; i < derivative.length - 1; i++) {
						if ((derivative[i] / derivative[i + 1] < 0) && derivative[i]>0){
//							double ratio = Math.abs((double) derivative[i]/ (double) derivative[i+1]);
//							System.out.println("RATIO "+i+" "+Math.abs(derivative[i]));
//							if ((threshold>0)&&(ratio<threshold))
							if ((threshold>0)&&(Math.abs(derivative[i])>threshold))
								d[i] = true;
						}
						else
							d[i] = false;
					}
					double max = getMax(derivative);
					if (max==derivative[derivative.length - 1])
						d[derivative.length - 1] = true;
					else
						d[derivative.length - 1] = false;
				}

				return d;
			}
			
		// returns a list of spikes indexes
		public static boolean[] findSpikes(double[] derivative,double threshold) {
			boolean[] d = new boolean[derivative.length];

			if (d.length > 0) {
				d[0] = false;
				for (int i = 1; i < derivative.length - 1; i++) {
					if (derivative[i] / derivative[i + 1] < 0){
//						double ratio = Math.abs((double) derivative[i]/ (double) derivative[i+1]);
//						System.out.println("RATIO "+i+" "+Math.abs(derivative[i]));
//						if ((threshold>0)&&(ratio<threshold))
						if ((threshold>0)&&(Math.abs(derivative[i])>threshold))
							d[i] = true;
					}
					else
						d[i] = false;
				}
				d[derivative.length - 1] = false;
			}

			return d;
		}

		// returns a list of spikes indexes
		public static boolean[] findSpikes(double[] derivative) {
			return findSpikes(derivative,-1);
		}
		
		

		// searches for an index into an array
		public static boolean isIn(List<Integer> indexarray, int index) {
			
			int size = indexarray.size();
			
			for (int i = 0; i < size; i++) {
				if (index == indexarray.get(i).intValue())
					return true;
			}
			
			return false;
		}
		
		
		// finds the indexes of zero points
		public static List<Integer> findZeros(double[] points) {
			
			int size = points.length;
			ArrayList<Integer> zeros = new ArrayList<Integer>();
			
			for (int i = 0; i < size; i++) {
				if (points[i]==0){
					int start = i;
					int end = i;
					
					for (int j=i+1;j<size;j++)
					{
						if (points[j]!=0){
							end = j-1;
							break;
						}
					}
					int center = start+((end-start)/2); 
					zeros.add(center);
					i = end;
				}
			}
			
			return zeros;
			
		}
		
		
		public static double[] logSubdivision(double start,double end,int numberOfParts){
			
			
			if (end<=start)
				return null;
			
			if (start == 0){
				start = 0.01;
			}
			double logStart = Math.log(start);
			double logEnd = Math.log(end);
			double step =0 ;
			if (numberOfParts >0){
				
				double difference = logEnd-logStart;
				step = (difference/(double)numberOfParts);
				
			}
//			double [] points = new double[numberOfParts+1];
			double[] linearpoints = new double[numberOfParts+1];
			
			for (int i=0;i<numberOfParts+1;i++){
				
//				points[i] = logStart+(i*step);
				
				linearpoints[i]= Math.exp(logStart+(i*step));
				if (linearpoints[i]<0.011)
					linearpoints[i] = 0;
			}
			
			return linearpoints;
		}
		
		
		public static double cohensKappaForDichotomy(long NumOf_A1_B1, long NumOf_A1_B0, long NumOf_A0_B1, long NumOf_A0_B0){
			long  T = NumOf_A1_B1+NumOf_A1_B0+NumOf_A0_B1+NumOf_A0_B0;
			
			double Pra = (double)(NumOf_A1_B1+NumOf_A0_B0)/(double) T ;
			double Pre1 = (double) (NumOf_A1_B1+NumOf_A1_B0) * (double) (NumOf_A1_B1+NumOf_A0_B1)/(double) (T*T);
			double Pre2 = (double) (NumOf_A0_B0+NumOf_A0_B1) * (double) (NumOf_A0_B0+NumOf_A1_B0)/(double) (T*T);
			double Pre = Pre1+Pre2;
			double Kappa = (Pra-Pre)/(1d-Pre);
			return roundDecimal(Kappa,3);
		}
		
		public static String kappaClassificationLandisKoch(double kappa){
			if (kappa<0)
				return "Poor";
			else if ((kappa>=0)&&(kappa<=0.20))
				return "Slight";
			else if ((kappa>=0.20)&&(kappa<=0.40))
				return "Fair";
			else if ((kappa>0.40)&&(kappa<=0.60))
				return "Moderate";
			else if ((kappa>0.60)&&(kappa<=0.80))
				return "Substantial";
			else if (kappa>0.80)
				return "Almost Perfect";
			else
				return "Not Applicable";
		}
		
		public static String kappaClassificationFleiss(double kappa){
			if (kappa<0)
				return "Poor";
			else if ((kappa>=0)&&(kappa<=0.40))
				return "Marginal";
			else if ((kappa>0.4)&&(kappa<=0.75))
				return "Good";
			else if (kappa>0.75)
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

		public static double getMax(double[] points) {
			double max = -Double.MAX_VALUE;
			for (int i = 0; i < points.length; i++) {
				if (max < points[i])
					max = points[i];
			}
			return max;
		}

		public static int getMax(int[] points) {
			int max = -Integer.MAX_VALUE;
			for (int i = 0; i < points.length; i++) {
				if (max < points[i])
					max = points[i];
			}
			return max;
		}

		public static int getMin(int[] points) {
			int min = Integer.MAX_VALUE;
			for (int i = 0; i < points.length; i++) {
				if (min > points[i])
					min = points[i];
			}
			return min;
		}

		public static double getMin(double[] points) {
			double min = Double.MAX_VALUE;
			for (int i = 0; i < points.length; i++) {
				if (min > points[i])
					min = points[i];
			}
			return min;
		}

		// calculates the frequency distribution for a set of points respect to a set of intervals
		public static double[] calcFrequencies(double[] interval, double[] points) {
			int intervs = interval.length;
			int npoints = points.length;
			double[] frequencies = new double[intervs];
			for (int i = 0; i < intervs; i++) {

				for (int j = 0; j < npoints; j++) {

					if (((i == 0) && (points[j] < interval[i])) || ((i == intervs - 1) && (points[j] >= interval[i - 1]) && (points[j] <= interval[i])) || ((i > 0) && (points[j] >= interval[i - 1]) && (points[j] < interval[i]))) {
						// System.out.println("(" + (i == 0 ? "" : interval[i - 1]) + "," + interval[i] + ")" + " - " + points[j]);
						frequencies[i] = frequencies[i] + 1;
					}
				}
			}

			return frequencies;
		}

		public static double[] normalizeFrequencies(double[] frequencies, int numberOfPoints) {
			int intervs = frequencies.length;
			for (int i = 0; i < intervs; i++) {
				frequencies[i] = frequencies[i] / (double) numberOfPoints;
			}

			return frequencies;

		}

		// checks if an interval contains at least one element from a sequence of points
		public static boolean intervalContainsPoints(double min, double max, double[] points) {
			// System.out.println(min+"-"+max);
			boolean contains = false;
			for (int i = 0; i < points.length; i++) {
				if ((points[i] >= min) && (points[i] < max)) {
					// System.out.println("---->"+points[i]);
					contains = true;
					break;
				}
			}
			return contains;
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
					if (!intervalContainsPoints(min + gap * i, rightmost, points)) {
						notcontains = true;
						break;
					}
				}

				// if there are empty intervals and there is space for another subdivision proceed
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

		public double[][] standardize(double[][] matrix) {
			return standardize(matrix, null, null);
		}

		public double[] means;
		public double[] SDs;

		// standardizes a matrix: each row represents a vector: outputs columns means and variances
		public double[][] standardize(double[][] matrix, double[] meansVec, double[] sdVec) {

				int ncols = matrix[0].length;
				int mrows = matrix.length;

				if ((means == null) && (SDs == null)) {
					means = new double[ncols];
					SDs = new double[ncols];
				}
				
				for (int i = 0; i< ncols; i++) {
					means[i] = columnMean(matrix, i);
					SDs[i] = columnSD(matrix, i);
				}
				
				double [][] matrix_std = new double[matrix.length][matrix[0].length];
				for (int i = 0; i< ncols; i++) {
					double mean_i=means[i];
					double sd_i=SDs[i];
					for (int j=0;j<mrows;j++) {
						
						matrix_std[j][i] = (matrix[j][i]-mean_i)/sd_i;
						
					}
				}
				
				return matrix_std;
			}

		// calculates the number of elements to take from a set with inverse weight respect to the number of elements
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

			double[] ret = { (-1d * b + Math.sqrt(b * b + 4 * a * (Math.abs(y) - c))) / (2 * a), (-1d * b - Math.sqrt(b * b + 4 * a * (Math.abs(y) - c))) / (2 * a) };
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

	
	
		//distributes uniformly elements in parts
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
			 * int numOfChunks = numberOfElements / chunksize; if ((numberOfElements % chunksize) != 0) numOfChunks += 1;
			 */

			return chunksize;
		}

		
		public static double[] uniformSampling(double min, double max, int maxElementsToTake){
			double step = (max-min)/(double)(maxElementsToTake-1);
			double [] samples = new double [maxElementsToTake];
			
			for (int i=0;i<samples.length;i++){
				double value = min+i*step;
				if (value>max)
					value=max;
				samples [i] = value;
			}
			
			return samples;
		}
		
		public static int[] uniformIntegerSampling(double min, double max, int maxElementsToTake){
			double step = (max-min)/(double)(maxElementsToTake-1);
			int [] samples = new int [maxElementsToTake];
			
			for (int i=0;i<samples.length;i++){
				double value = min+i*step;
				if (value>max)
					value=max;
				samples [i] = (int)value;
			}
			
			return samples;
		}
		
	
		public static ExampleSet matrix2ExampleSet(double[][] sampleVectors) {

			int m = sampleVectors.length;

			BigSamplesTable samples = new BigSamplesTable();

			for (int k = 0; k < m; k++)
				samples.addSampleRow("sample", sampleVectors[k]);

			return samples.generateExampleSet();

		}

		public static double[][] exampleSet2Matrix(ExampleSet set) {

			int m = set.size();
			ExampleTable table = set.getExampleTable();
			int n = table.getAttributeCount();
			double[][] matrix = new double[m][n - 1];
			for (int i = 0; i < m; i++) {
				DataRow row = table.getDataRow(i);
				for (int j = 0; j < n - 1; j++) {
					if (!table.getAttribute(j).isNominal()) {
						double d = row.get(table.getAttribute(j));
						matrix[i][j] = d;
					}
				}
			}

			return matrix;

		}

		// gets all the columns from a matrix
		public static double[][] traspose(double[][] matrix) {
			int m = matrix.length;
			if (m > 0) {
				int n = matrix[0].length;

				double columns[][] = new double[n][m];

				for (int i = 0; i < n; i++) {
					for (int j = 0; j < m; j++)
						columns[i][j] = matrix[j][i];
				}

				return columns;
			} else
				return null;
		}

		// gets a column from a matrix
		public static double[] getColumn(int index, double[][] matrix) {
			int colulen = matrix.length;
			double column[] = new double[colulen];
			for (int i = 0; i < colulen; i++) {
				column[i] = matrix[i][index];
			}
			return column;
		}

		// gets a column from a matrix
		public static void substColumn(double[] column, int index, double[][] matrix) {

			for (int i = 0; i < matrix.length; i++) {
				matrix[i][index] = column[i];
			}

		}

		// merge matrixes: puts the rows of a matrix under another matrix
		public static double[][] mergeMatrixes(double[][] matrix1, double[][] matrix2) {

			if ((matrix1 == null) || (matrix1.length == 0))
				return matrix2;
			else if ((matrix2 == null) || (matrix2.length == 0))
				return matrix1;
			else {
				int len1 = matrix1.length;
				int len2 = matrix2.length;
				int superlen = len1 + len2;
				double[][] supermatrix = new double[superlen][];
				for (int i = 0; i < len1; i++) {
					supermatrix[i] = matrix1[i];
				}
				for (int i = len1; i < superlen; i++) {
					supermatrix[i] = matrix2[i - len1];
				}
				return supermatrix;
			}
		}

		public static String vector2String(double[] vector) {
			String out = "";
			for (int i = 0; i < vector.length; i++) {
				if (i > 0)
					out = out + "," + vector[i];
				else
					out = "" + vector[i];
			}

			return out;
		}


		public static double indexString(String string) {
			// string = Sha1.SHA1(string);
			StringBuffer sb = new StringBuffer();
			if ((string == null) || (string.length() == 0))
				return -1;

			int m = string.length();
			for (int i = 0; i < m; i++) {
				sb.append((int) string.charAt(i));
			}

			double d = Double.MAX_VALUE;
			try {
				d = Double.valueOf(sb.toString());
			} catch (Throwable e) {
			}

			if (d > Integer.MAX_VALUE)
				return (indexString(string.substring(0, 3)));

			return d;
		}
		
		
		public static double [] featureTimesInSec(double windowShift, File audio) throws Exception{
			AudioBits bits = new AudioBits(audio);
			short[] signal = bits.getShortVectorAudio();
			bits.ais.close();
			float sfrequency = bits.getAudioFormat().getSampleRate();
			double total_audio_sec = (double)signal.length/(double) sfrequency; 
			double t = 0;
			int ntimes = (int)Math.floor(total_audio_sec/windowShift);
			double times [] = new double[ntimes];
			int counter = 0;
			while (t<total_audio_sec) {
				if (counter<ntimes)
					times[counter] = t;
				t = t+windowShift;
				counter++;
			}
			return times;
			
		}

		public static double [] featureTimesInSec(double windowShift, int sfrequency, int signalLength) throws Exception{
			double total_audio_sec = (double)signalLength/(double) sfrequency; 
			double t = 0;
			int ntimes = (int)Math.floor(total_audio_sec/windowShift);
			double times [] = new double[ntimes];
			int counter = 0;
			while (t<total_audio_sec) {
				if (counter<ntimes)
					times[counter] = t;
				t = t+windowShift;
				counter++;
			}
			return times;
			
		}
		
		/*
		public static List<String> parseCVSString(String row, String delimiter) throws Exception {

			List<String> elements = new ArrayList<String>();
			String phrase = row;
			int idxdelim = -1;
			boolean quot = false;
			phrase = phrase.trim();
			while ((idxdelim = phrase.indexOf(delimiter)) >= 0) {
				quot = phrase.startsWith("\"");
				if (quot) {
					phrase = phrase.substring(1);
					String quoted = "";
					if (phrase.startsWith("\""))
						phrase = phrase.substring(1);
					else{
						RE regexp = new RE("[^\\\\]\"");
						boolean matching = regexp.match(phrase);

						if (matching) {
							int i0 = regexp.getParenStart(0);
							quoted = phrase.substring(0, i0 + 1).trim();
							phrase = phrase.substring(i0 + 2).trim();
						}
					}

					if (phrase.startsWith(delimiter))
						phrase = phrase.substring(1);

					elements.add(quoted);

				} else {
					elements.add(phrase.substring(0, idxdelim));
					phrase = phrase.substring(idxdelim + 1).trim();
				}
			}
			if (phrase.startsWith("\""))
				phrase = phrase.substring(1);

			if (phrase.endsWith("\""))
				phrase = phrase.substring(0, phrase.length() - 1);

			elements.add(phrase);

			return elements;
		}
		*/
		
		public static double[] quantiles (double[] data) {
			
		        Arrays.sort(data);
		        double q1 = calculateQuartile(data, 0.25);
		        double q2 = calculateQuartile(data, 0.5);
		        double q3 = calculateQuartile(data, 0.75);
		        
		        double [] quants = new double[3];
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
		 
		 
		 public static double[] initializeVector(int n, double value) {
		        double[] numericArray = new double[n];
		        for (int i = 0; i < n; i++) {
		            numericArray[i] = value;
		        }

		        return numericArray;
		    }
		 
		 public static double[] reduceVectorValues(double values[] , double reductionFactor) {
		        double[] numericArray = new double[values.length];
		        
		        for (int i = 0; i < values.length; i++) {
		            numericArray[i] = values[i]-reductionFactor*values[i];
		        }

		        return numericArray;
		    }
		
}
