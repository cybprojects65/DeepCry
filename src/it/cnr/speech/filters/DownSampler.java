package it.cnr.speech.filters;

public class DownSampler {

	
	    public static double[] downsample(double[] input, int factor) {
	    	
	        int outputLength = input.length / factor;
	        
	        double[] output = new double[outputLength];

	        for (int i = 0; i < outputLength; i++) {
	            output[i] = input[i * factor];
	        }

	        //System.out.println("Original Signal: " + Arrays.toString(input));
	        //System.out.println("Downsampled Signal: " + Arrays.toString(output));
	        
	        return output;
	    }

	    
	
	
}
