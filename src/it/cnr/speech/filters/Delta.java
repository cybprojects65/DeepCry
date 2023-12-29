package it.cnr.speech.filters;

public class Delta {

	public static void calcDelta(double A[][], int numCoeff) throws Exception {
		int delta;

		for (int j = 0; j < numCoeff; j++) {
			delta = j + numCoeff;
			completeDelta(A, j, delta);
		}
	}

	public static void calcDoubleDelta(double A[][], int numCoeff) throws Exception {
		int fine = numCoeff * 2;
		for (int delta = numCoeff; delta < fine; delta++) {
			int doppioDelta = delta + numCoeff;
			completeDelta(A, delta, doppioDelta);
		}
	}

	private static void completeDelta(double A[][], int j, int d) throws Exception {
		if (A.length < 4) {
			throw new Exception();
		}
		if (A.length > 2) {
			A[0][d] = A[1][j] - A[0][j];
			A[1][d] = (A[2][j] - A[0][j]) - (A[0][j] - A[1][j]) / 4;
		}
		for (int i = 2; i < A.length - 2; i++)
			A[i][d] = ((2 * (A[i + 2][j] - A[i - 2][j])) + (A[i + 1][j] - A[i - 1][j])) / 8;
		if (A.length > 3) {
			A[A.length - 2][d] = (A[A.length - 1][j] - A[A.length - 3][j]) - (A[A.length - 1][j] - A[A.length - 2][j]);
			A[A.length - 1][d] = A[A.length - 2][j] - A[A.length - 1][j];
		}
	}

}
