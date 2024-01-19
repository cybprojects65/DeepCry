package it.cnr.clustering;

import java.util.LinkedList;
import java.util.List;

import com.rapidminer.example.Attribute;
import com.rapidminer.example.ExampleSet;
import com.rapidminer.example.table.AttributeFactory;
import com.rapidminer.example.table.DoubleArrayDataRow;
import com.rapidminer.example.table.MemoryExampleTable;
import com.rapidminer.tools.Ontology;

public abstract class SamplesTable {

	// attributes = columns of numbers
	public double minY = 0;
	public double maxY = 0;
	public double minX = 0;
	public double maxX = 0;

	public ExampleSet generateExampleSet() {
		// create attribute list
		List<Attribute> attributes = new LinkedList<Attribute>();
		// generate columns for attributes
		for (int a = 0; a < getNumOfAttributes(); a++) {
			attributes.add(AttributeFactory.createAttribute("att" + a, Ontology.REAL));
		}

		// add a label column
		Attribute label = AttributeFactory.createAttribute("label", Ontology.NOMINAL);
		attributes.add(label);

		// create table
		MemoryExampleTable table = new MemoryExampleTable(attributes);

		// fill table (here : only real values )
		for (int d = 0; d < getNumOfDataRows(); d++) {
			// generate a row of double values
			double[] data = new double[attributes.size()];
			// fill rows data
			for (int a = 0; a < getNumOfAttributes(); a++) {
				// all with proper data here
				data[a] = getValue(d, a);
			}
			// maps the nominal classifcation to a double value
			data[data.length - 1] = label.getMapping().mapString(getClassification(d));

			// add data row
			table.addDataRow(new DoubleArrayDataRow(data));
		}

		// create example set
		ExampleSet exampleSet = table.createExampleSet(label);

		return exampleSet;

	}

	public void generateSampleTable(ExampleSet es) {

		MemoryExampleTable met = (MemoryExampleTable) es.getExampleTable();
		int numofcolumns = met.getAttributeCount();
		int numofrows = met.size();
		// System.out.println("COL "+numofcolumns+" ROWS "+numofrows);

		for (int i = 0; i < numofrows; i++) {

			Attribute labelAtt = met.getAttribute(numofcolumns - 1);
			int index = (int) met.getDataRow(i).get(labelAtt);
			String label = labelAtt.getMapping().mapIndex(index);
			addSampleRow(label, 0);

			// addLabel(i,label);

			for (int j = 0; j < numofcolumns - 1; j++) {
				Attribute att = AttributeFactory.createAttribute("att" + j, Ontology.REAL);
				att.setTableIndex(j);
				// System.out.println("ADDING TO " + i+","+j);
				DoubleArrayDataRow dadr = (DoubleArrayDataRow) met.getDataRow(i);
				double element = dadr.get(att);

				addSample(i, j, element);
			}

		}

	}
	
	
	abstract public int getNumOfAttributes();

	abstract public int getNumOfDataRows();

	abstract public double getValue(int d, int a);

	abstract public String getClassification(int d);

	public String toString() {

		StringBuffer bs = new StringBuffer();

		bs.append("NUMBER OF ROWS: " + getNumOfDataRows() + "\n");
		bs.append("NUMBER OF COLUMNS: " + getNumOfAttributes() + "\n");

		for (int i = 0; i < getNumOfDataRows(); i++) {
			bs.append("ROW " + i + " : ");
			bs.append("LABEL " + getClassification(i) + " : ");

			for (int j = 0; j < getNumOfAttributes(); j++) {
				bs.append(getValue(i, j) + "\t");
			}

			bs.append("\n");
		}

		return bs.toString();
	}

	public void calculateBounds() {

		int Ylen = getNumOfAttributes();
		int Xlen = getNumOfDataRows();

		for (int i = 0; i < Xlen; i++) {
			for (int j = 0; j < Ylen; j++) {
				double localmin = minY;
				double localmax = maxY;
				if (j == 0) {
					localmin = minX;
					localmax = maxX;
				}
				double point = getValue(i, j);
				if (point < localmin) {
					localmin = point;
				} else if (point > localmax) {
					localmax = point;
				}
				if (j == 0) {
					minX = localmin;
					maxX = localmax;
				} else {
					minY = localmin;
					maxY = localmax;
				}
			}
		}

	}
	
	abstract public void addLabel(int i, String label);

	abstract public void addSample(int i, int j, double value);

	abstract public void addSampleRow(String label, double... values);

}

