package it.cnr.evaluation;

public class Range {
	
	public double t0;
	public double t1;
	
	public boolean isInRange(double t) {
		
		return (t>=t0 && t<=t1);
		
	}
	
	public Range(double t0, double t1) {
		this.t0 = t0;
		this.t1 = t1;
	}

	public String toString() {
		return t0+";"+t1;
	}
}
