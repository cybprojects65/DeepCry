package it.cnr.speech.filters;

public class HalfWaveRectification {
    public static double[] applyHalfWaveRectification(double[] input) {
        int N = input.length;

        double[] output = new double[N];
        for (int i = 0; i < N; i++) {
                output[i] = Math.max(0, input[i]);
        }

        return output;
    }

    
}