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

public class UtilsObjects {
	
	// ###########Java object manipulation
		public static void saveObject(File outputFile, Object o) throws Exception {
			ObjectOutputStream oos = new ObjectOutputStream(new FileOutputStream(outputFile));
			oos.writeObject(o);
			oos.close();
		}

		public static Object loadObject(File inputFile) throws Exception {
			ObjectInputStream ois = new ObjectInputStream(new FileInputStream(inputFile));
			return (ois.readObject());
		}

		// ###########Data preparation
		public static int saveTS4DL(List<double[][]> timeSeries, int label, int counter, File trainingFeaturesFolder,
				File testFeaturesFolder, File trainingLabelsFolder, File testLabelsFolder) throws Exception {

			for (double[][] ts : timeSeries) {
				StringBuffer sb = new StringBuffer();
				for (int i = 0; i < ts.length; i++) {
					double[] row = ts[i];
					for (int j = 0; j < row.length; j++) {
						double r = row[j];
						sb.append(r);
						if (j < row.length - 1)
							sb.append(",");
					}
					sb.append("\n");
				}

				File featureFileTr = new File(trainingFeaturesFolder, counter + ".csv");
				File featureFileTest = new File(testFeaturesFolder, counter + ".csv");
				FileWriter fwTr = new FileWriter(featureFileTr);
				fwTr.write(sb.toString().trim());
				fwTr.close();
				FileWriter fwTest = new FileWriter(featureFileTest);
				fwTest.write(sb.toString().trim());
				fwTest.close();

				File labelFileTr = new File(trainingLabelsFolder, counter + ".csv");
				File labelFileTest = new File(testLabelsFolder, counter + ".csv");
				FileWriter lTr = new FileWriter(labelFileTr);
				lTr.write("" + label);
				lTr.close();
				FileWriter lTest = new FileWriter(labelFileTest);
				lTest.write("" + label);
				lTest.close();

				counter++;
			}

			return counter;

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

}
