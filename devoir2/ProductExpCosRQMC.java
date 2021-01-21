package devoir2;

import java.io.*;
import umontreal.ssj.hups.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;

// Question 3 du devoir 2

public class ProductExpCosRQMC implements MonteCarloModelDouble {

	int s;
	double b;
	double prod;

	// Constructor.
	public ProductExpCosRQMC(int s, double b) {
		this.s = s; 	this.b = b;
	}

	// Generates and returns X, without IS.
	public void simulate (RandomStream stream) {
		prod = 0.0;
		double u;
		for (int j = 0; j < s; j++) {
			u = stream.nextDouble();
			prod += u * Math.cos(b * u);
		}
		
		prod = Math.pow(prod,2);
	}

	// Generates and returns X, without IS.
	public double getPerformance () {
		return prod;
	}

	// Descriptor of this model.
	public String toString () {
		return "Test function for MC and RQMC: product of exponentials and cosine functions.";
	}

	public static void main(String[] args) throws IOException {
		int s = 3;
		int n = (int)Math.pow(2, 16);
		int m = 20;                     // Number of RQMC randomizations.
		RandomStream stream = new LFSR113();
		DigitalNet p = new SobolSequence(16, 31, s); // n = 2^{16} points in s dim.
		// PointSetRandomization rand = new LMScrambleShift(stream);
		PointSetRandomization rand = new RandomShift(stream);

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(new ProductExpCosRQMC(s, 0.5), stream, p, rand, n, m));
		
		System.out.print("\n\n\n");

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(new ProductExpCosRQMC(s, 50), stream, p, rand, n, m));
		
	}
}

