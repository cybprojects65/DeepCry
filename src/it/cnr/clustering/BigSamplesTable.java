package it.cnr.clustering;

import java.util.HashMap;

public class BigSamplesTable extends SamplesTable {

	BigSparseTable table;
	HashMap<Integer, String> classifications;
	Integer currentIndex;

	public BigSamplesTable() {
		table = new BigSparseTable();
		classifications = new HashMap<Integer, String>();
		currentIndex = 0;
	}

	@Override
	public int getNumOfAttributes() {
		return table.width().intValue();
	}

	@Override
	public int getNumOfDataRows() {

		return table.size().intValue();

	}

	@Override
	public double getValue(int d, int a) {

		return table.get(d, a);

	}

	@Override
	public String getClassification(int d) {

		return classifications.get(d);
	}

	public void addSampleRow(String label, double... values) {
		if (values != null) {
		classifications.put(currentIndex, label);
		int j = 0;
		for (Double value : values) {

			table.add(currentIndex, j, value);

			j++;
		}

		currentIndex = currentIndex + 1;
		}
	}

	public void addSample(int i, int j, double value) {

		if (i < currentIndex)
			table.add(i, j, value);
	}

	public void addLabel(int i, String label) {

		if (i < currentIndex)
			classifications.put(i, label);
	}

	
	public static void main(String[] args){
		
			BigSamplesTable bst = new BigSamplesTable();
			bst.addSampleRow("prova 1", 10, 12,13,14,15);
			bst.addSampleRow("prova 2", 20, 15,14,15);
			bst.addSampleRow("prova 3", 30, 11,110,150);
			bst.addSample(0, -1,150);
			System.out.println(bst.toString());
			
			bst.generateExampleSet();
			
			
	}
	
}
