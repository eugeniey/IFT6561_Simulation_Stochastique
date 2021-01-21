package ProbaRuine;

import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.*;
import umontreal.ssj.util.*;

public class ruineIS {

static final double lambda     = 1.0;    // Arrival rate of claims.
static final double alpha      = 0.5;    // 1 / expected claim size.
static final double beta   = 0.25;    // Arrival rate of claims.
static final double r0         = 200.0;  // Initial reserve.
static final double RE         = 0.01;  // Relative error.
double theta;

RandomVariateGen             genArrivals;   // For claim arrivals.
GammaRejectionLoglogisticGen genSizes;      // For claim sizes.
Tally statIs = new Tally ("Ruin probability with IS");

public double simulRuin (double c) {
   double sum = 0.0;
   while (sum < r0)
      sum += genSizes.nextDouble() - c * genArrivals.nextDouble();
   return Math.exp (- theta * sum); 
}

public ruineIS (double c, int n) {
   // Computes IS parameters and makes n simulation runs with IS,
   // with input rate c for the premiums.
   theta = (c - 8.0 + Math.sqrt(16.0 * c + c * c)) / (8.0 * c);
   double lambdaIs = lambda + c * theta; 
   double betaIs = beta - theta;
   genArrivals = new RandomVariateGen
      (new MRG32k3a(), new ExponentialDist (lambdaIs));
   genSizes = new GammaRejectionLoglogisticGen
      (new MRG32k3a(),new MRG32k3a(), alpha, betaIs);
   for (int i=0; i < n; i++) {
      statIs.add (simulRuin (c));
      }
   System.out.println (" lambda = " + 
      PrintfFormat.format (8, 3, 1, lambda));
   System.out.println (" lambdaIS = " + 
		   PrintfFormat.format (8, 3, 1, lambdaIs));
   System.out.println (" alpha   = " + 
      PrintfFormat.format (8, 3, 1, alpha));
   System.out.println (" beta   = " + 
		   PrintfFormat.format (8, 3, 1, beta));
   System.out.println (" betaIS   = " + 
		   PrintfFormat.format (8, 3, 1, betaIs));
   System.out.println("theta = " + theta);
   System.out.println (" c      = " + 
      PrintfFormat.format (8, 3, 1, c));
   System.out.println (" R(0)   = " + 
      PrintfFormat.format (8, 3, 1, r0));
   System.out.println (" n      = " + n);
   System.out.println ();
   double var = statIs.variance();
   //double sigma = statIs.standardDeviation();
   System.out.println (statIs.formatCIStudent (0.90));
   System.out.println (" Variance with IS = " + PrintfFormat.format(10, 2, 2, var));
   double p = statIs.average();
   System.out.println(" p      = " + PrintfFormat.format(10,2,2,p));
   //System.out.println("Relative error   = " + PrintfFormat.format(10, 5, 2, (sigma/p)));
   System.out.println (" Sample size for 1% error with MC = " + 
       PrintfFormat.format( 10, 2, 2, 10000.0 * (1.0-p) / p));
   System.out.println (" Sample size for 1% error with IS = " + 
      PrintfFormat.format( 10, 2, 2, 10000.0 * statIs.variance() / (p*p)));
   System.out.println ();
   System.out.println ("---------------------------------------------");
}

public static void main (String[] args) { 
   new ruineIS (3.0, 10000);
   new ruineIS (5.0, 10000);
}
}
