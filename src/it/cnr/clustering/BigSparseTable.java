package it.cnr.clustering;

import java.util.HashMap;

public class BigSparseTable {

	HashMap<Integer, HashMap<Integer, Double>> table;
	Integer tableSize;
	Integer tableWidth;

	public BigSparseTable() {
		table = new HashMap<Integer, HashMap<Integer, Double>>();
		tableSize = 0;
		tableWidth = 0;
	}

	public Integer size() {

		return tableSize;

	}

	public Integer width() {

		return tableWidth;

	}

	public void add(Integer i, Integer j, double value) {

//		System.out.println("ADDING " + i + "," + j);
		if ((i<0)||(j<0))
			return;
		
		double val = get(i, j);

		if (val != 0) {
			table.get(i).put(j, value);
		}

		else {
			HashMap<Integer, Double> row;

			// if size<=i create a new hashmap
			if (tableSize <= i) {
				row = new HashMap<Integer, Double>();
				table.put(i, row);

				tableSize = i + 1;
			}
			// else get i-th hashmap
			else {
				row = table.get(i);
				if (row == null) {
					row = new HashMap<Integer, Double>();
					table.put(i, row);
				}
			}

			row.put(j, value);

			if (tableWidth <= j)
				tableWidth = j + 1;
		}
	}

	// default is 0
	public double get(Integer i, Integer j) {

		Double value = null;

		if (tableSize.compareTo(i) > 0) {
			value = table.get(i).get(j);
		}

		if (value == null)
			value = Double.valueOf(0);

		return value;
	}

}
