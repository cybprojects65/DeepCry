package it.cnr.workflow.utils;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.ObjectInputStream;
import java.io.ObjectOutputStream;
import java.util.ArrayList;
import java.util.List;

import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.DataRow;
import com.rapidminer.example.table.ExampleTable;

import it.cnr.clustering.BigSamplesTable;
import it.cnr.speech.audiofeatures.AudioBits;

public class UtilsVectorMatrix {

	// ###########Matrix/Vector manipulation
	public static double[][] subsetRows(double[][] matrix, int row0, int row1) {

		double[][] sub = new double[row1 - row0 + 1][matrix.length];

		for (int i = row0; i <= row1; i++) {

			sub[i - row0] = matrix[i];

		}

		return sub;

	}

	public static double[] columnToVector(double[][] matrix, int i) {

		double[] vector = new double[matrix.length];
		for (int k = 0; k < matrix.length; k++) {
			vector[k] = matrix[k][i];
		}
		return vector;
	}

	public static double columnMean(double[][] matrix, int i) {

		double[] vector = columnToVector(matrix, i);
		return UtilsMath.mean(vector);
	}

	public static double columnSD(double[][] matrix, int i) {

		double[] vector = columnToVector(matrix, i);
		double variance = com.rapidminer.tools.math.MathFunctions.variance(vector, Double.NEGATIVE_INFINITY);
		return Math.sqrt(variance);
	}

	public static double[] columnQ1(double[][] matrix) {

		int ncol = matrix[0].length;
		double[] o = new double[ncol];

		for (int j = 0; j < ncol; j++) {
			double[] vector = columnToVector(matrix, j);
			double q1 = UtilsMath.quantiles(vector)[0];
			o[j] = q1;
		}
		return o;
	}

	public static double[] columnQ3(double[][] matrix) {

		int ncol = matrix[0].length;
		double[] o = new double[ncol];

		for (int j = 0; j < ncol; j++) {
			double[] vector = columnToVector(matrix, j);
			double q3 = UtilsMath.quantiles(vector)[2];
			o[j] = q3;
		}
		return o;
	}

	public static double[] columnMeans(double[][] matrix) {

		int ncol = matrix[0].length;
		double[] colmeans = new double[ncol];

		for (int j = 0; j < ncol; j++) {

			double mean = columnMean(matrix, j);
			colmeans[j] = mean;
		}
		return colmeans;
	}

	public static double[] columnSDs(double[][] matrix) {

		int ncol = matrix[0].length;
		double[] colsds = new double[ncol];

		for (int j = 0; j < ncol; j++) {
			double sd = columnSD(matrix, j);
			colsds[j] = sd;
		}

		return colsds;
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
				if (points[i] == 0) {
					int start = i;
					int end = i;

					for (int j = i + 1; j < size; j++) {
						if (points[j] != 0) {
							end = j - 1;
							break;
						}
					}
					int center = start + ((end - start) / 2);
					zeros.add(center);
					i = end;
				}
			}

			return zeros;

		}
	public static double getArgMax(double[] points) {
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

	public static double getArgMin(double[] points) {
			double min = Double.MAX_VALUE;
			for (int i = 0; i < points.length; i++) {
				if (min > points[i])
					min = points[i];
			}
			return min;
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

	public double[] means;
	public double[] SDs;

	// standardizes a matrix: each row represents a vector: outputs columns means
	// and variances
	public double[][] standardize(double[][] matrix, double[] meansVec, double[] sdVec) {

		int ncols = matrix[0].length;
		int mrows = matrix.length;

		if ((means == null) && (SDs == null)) {
			means = new double[ncols];
			SDs = new double[ncols];
		}

		for (int i = 0; i < ncols; i++) {
			means[i] = columnMean(matrix, i);
			SDs[i] = columnSD(matrix, i);
		}

		double[][] matrix_std = new double[matrix.length][matrix[0].length];
		for (int i = 0; i < ncols; i++) {
			double mean_i = means[i];
			double sd_i = SDs[i];
			for (int j = 0; j < mrows; j++) {

				matrix_std[j][i] = (matrix[j][i] - mean_i) / sd_i;

			}
		}

		return matrix_std;
	}

	public double[][] standardize(double[][] matrix) {
		return standardize(matrix, null, null);
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

	// substitutes a column in a matrix
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
		
	public static double[] initializeVector(int n, double value) {
		double[] numericArray = new double[n];
		for (int i = 0; i < n; i++) {
			numericArray[i] = value;
		}

		return numericArray;
	}

	public static double[] reduceVectorValues(double values[], double reductionFactor) {
		double[] numericArray = new double[values.length];

		for (int i = 0; i < values.length; i++) {
			numericArray[i] = values[i] - reductionFactor * values[i];
		}

		return numericArray;
	}
	
	

	/*
	 * public static List<String> parseCVSString(String row, String delimiter)
	 * throws Exception {
	 * 
	 * List<String> elements = new ArrayList<String>(); String phrase = row; int
	 * idxdelim = -1; boolean quot = false; phrase = phrase.trim(); while ((idxdelim
	 * = phrase.indexOf(delimiter)) >= 0) { quot = phrase.startsWith("\""); if
	 * (quot) { phrase = phrase.substring(1); String quoted = ""; if
	 * (phrase.startsWith("\"")) phrase = phrase.substring(1); else{ RE regexp = new
	 * RE("[^\\\\]\""); boolean matching = regexp.match(phrase);
	 * 
	 * if (matching) { int i0 = regexp.getParenStart(0); quoted =
	 * phrase.substring(0, i0 + 1).trim(); phrase = phrase.substring(i0 + 2).trim();
	 * } }
	 * 
	 * if (phrase.startsWith(delimiter)) phrase = phrase.substring(1);
	 * 
	 * elements.add(quoted);
	 * 
	 * } else { elements.add(phrase.substring(0, idxdelim)); phrase =
	 * phrase.substring(idxdelim + 1).trim(); } } if (phrase.startsWith("\""))
	 * phrase = phrase.substring(1);
	 * 
	 * if (phrase.endsWith("\"")) phrase = phrase.substring(0, phrase.length() - 1);
	 * 
	 * elements.add(phrase);
	 * 
	 * return elements; }
	 */

	//binarises a matrix based on thresholds per column
	public static double[][] binariseMatrix(double [][] matrix, double[] binarisationThresholds){
		int nrow = matrix.length;
		int ncol = matrix[0].length;
		double [][] binarymatrix = new double[nrow][ncol];
		
		for (int i=0;i<nrow;i++) {
			for (int j=0;j<ncol;j++) {
				if (matrix[i][j]>binarisationThresholds[j]) {
					binarymatrix[i][j] = 1;
				}else
					binarymatrix[i][j] = 0;
			}
		}
		
		return binarymatrix;
	}

	

}
