package it.cnr.clustering;

import java.io.File;
import java.io.FileWriter;
import java.util.HashMap;
import java.util.LinkedHashMap;

import com.rapidminer.RapidMiner;
import com.rapidminer.example.Attribute;
import com.rapidminer.example.Attributes;
import com.rapidminer.example.Example;
import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.DataRow;
import com.rapidminer.example.table.ExampleTable;
import com.rapidminer.operator.IOContainer;
import com.rapidminer.operator.IOObject;
import com.rapidminer.operator.clustering.Cluster;
import com.rapidminer.operator.clustering.ClusterModel;
import com.rapidminer.tools.OperatorService;

public class KMeans {

	int minPoints;

	HashMap<String, String> labels;
	public LinkedHashMap<String, Integer> pointsPerCluster = new LinkedHashMap<>();
	public static Boolean RapidMinerInitialised = false;

	public void initRapidMiner() {
		if (!RapidMinerInitialised) {
			System.setProperty("rapidminer.init.operators", "./config/operators.xml");
			RapidMiner.init();
			System.out.println("Rapid Miner initialized");
			RapidMinerInitialised = true;
		}
	}

	public KMeans(HashMap<String, String> labels) {

		this.labels = labels;
		initRapidMiner();
	}

	public File compute(int kk, int maxRuns, int maxOptimizations, int minPoints, double[][] features, File outputFile)
			throws Exception {
		this.minPoints = minPoints;

		// take elements and produce example set
		com.rapidminer.operator.clustering.clusterer.KMeans kmeans = (com.rapidminer.operator.clustering.clusterer.KMeans) OperatorService
				.createOperator("KMeans");

		kmeans.setParameter("k", "" + kk);
		kmeans.setParameter("max_runs", "" + maxRuns);
		kmeans.setParameter("max_optimization_steps", "" + maxOptimizations);
		kmeans.setParameter("keep_example_set", "true");
		kmeans.setParameter("add_cluster_attribute", "true");

		ExampleSet points = matrix2ExampleSet(features);
		IOContainer innerInput = new IOContainer(points);

		//System.out.println("KMeans: Clustering...");
		long ti = System.currentTimeMillis();
		IOContainer output = kmeans.apply(innerInput);
		//System.out.println("KMEANS: ...ELAPSED CLUSTERING TIME: " + (System.currentTimeMillis() - ti));
		//System.out.println("KMeans: ...Clustering Finished");
		IOObject[] outputvector = output.getIOObjects();
		File outputTable = buildClusterTable(outputvector, outputFile);
		return outputTable;
	}

	protected File buildClusterTable(IOObject[] outputvector, File outputFile) throws Exception {

		ClusterModel innermodel = (ClusterModel) outputvector[0];
		ExampleSet es = (ExampleSet) outputvector[1];

		// System.out.println("Analyzing Cluster ->"+" minpoints"+minPoints);
		int nClusters = innermodel.getClusters().size();
		//System.out.println("Analyzing Cluster -> K = " + nClusters + " clusters");
		//System.out.println("Start Write On DB");
		StringBuffer bufferRows = new StringBuffer();
		int allCounter = 0;
		bufferRows.append("id,label,cluster_id,is_an_outlier\n");

		for (Cluster c : innermodel.getClusters()) {
			// take cluster id
			int id = c.getClusterId();
			boolean outlier = false;
			// take cluster element indexes
			int npoints = c.getExampleIds().size();
			// System.out.println("Analyzing Cluster ->"+id+" with "+npoints);
			pointsPerCluster.put("" + id, npoints);
			if (npoints <= minPoints)
				outlier = true;

			int k = 0;
			for (Object o : c.getExampleIds()) {
				// transform into a numerical index
				int idd = (int) Double.parseDouble("" + o);

				// take the corresponding sample
				Example e = es.getExample(idd - 1);
				// take the attributes of the sample
				Attributes attributes = e.getAttributes();

				// for each attribute (yet filtered on numeric ones) add to the writing row

				StringBuffer valueStrings = new StringBuffer();
				for (Attribute attribute : attributes) {
					valueStrings.append(e.getValue(attribute) + ",");
				}
				String towrite = valueStrings.toString();
				towrite = towrite.substring(0, towrite.length() - 1);

				// append the clusterid and outlier
				bufferRows.append(idd + "," + labels.get("" + idd) + "," + id + "," + outlier + "\n");

				k++;
				allCounter++;
//				logger.trace("DBScan: Classification : "+towrite+"->"+id+" is outlier?"+outlier);
			}
		}
		if (bufferRows.length() > 0) {

			//System.out.println("Writing into DB");
			File output = outputFile;
			FileWriter fw = new FileWriter(output);
			fw.write(bufferRows.toString());
			fw.close();

			//System.out.println("Finished with writing into DB");
			return output;
		} else {
			System.out.println("Nothing to write in the buffer");
			return null;
		}

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
