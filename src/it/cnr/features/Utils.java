package it.cnr.features;

public class Utils {

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
	
	public static double timeToSamples(double time, double fs) {
		
		return Math.round(fs * time);
		
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


}
