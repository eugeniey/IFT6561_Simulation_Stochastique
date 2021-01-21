package travailFinal;

import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.mcqmctools.*;
import umontreal.ssj.hups.BakerTransformedPointSet;
import umontreal.ssj.hups.DigitalNet;
import umontreal.ssj.hups.KorobovLattice;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.PointSetRandomization;
import umontreal.ssj.hups.RandomShift;
import umontreal.ssj.hups.SobolSequence;

import umontreal.ssj.util.RootFinder;

import umontreal.ssj.functions.*;

import java.lang.Math;

import devoir2.ProductExpCosRQMC;
import travailFinal.AsianPriceProcess;

/*
 * Inspired from class exemples
 */
public class AsianOption implements MonteCarloModelDouble {

	StochasticProcess priceProcess; // Underlying process for the price.
	int d; // Number of observation times.
	double[] obsTimes; // obsTimes[0..d] must contain obsTimes[0]=0.0,
	                   // plus the d positive observation times.	
	double[] path; // Sample path of the process.
	double strike; // Strike price.
	double discount; // Discount factor exp(-r * obsTimes[t]).


	public AsianOption(double r, int d, double[] obsTimes, double strike) {
		this.d = d;
		this.obsTimes = new double[d + 1];
		for (int j = 0; j <= d; j++)
			this.obsTimes[j] = obsTimes[j];
		this.strike = strike;
		discount = Math.exp(-r * obsTimes[d]);
	}

    // This constructor also specifies the underlying process.
	public AsianOption(StochasticProcess sp, double r, int d,
			double[] obsTimes, double strike) {
		this(r, d, obsTimes, strike);
		setProcess(sp);
	}


	public AsianOption(double r, int d, double T1, double T, double strike) {
		this.d = d;
		obsTimes = new double[d + 1];
		obsTimes[0] = 0.0;
		for (int j = 1; j <= d; j++)
			obsTimes[j] = T1 + (double) (j - 1) * (T - T1) / (double) (d - 1);
		this.strike = strike;
		discount = Math.exp(-r * obsTimes[d]);
	}

	public void setProcess(StochasticProcess sp) {
		// Reset the process to sp. Assumes that obsTimes have been set.
		priceProcess = sp;
		sp.setObservationTimes(obsTimes, d);
	}

	public double getPerformance() {
		double average = 0.0; // Average over sample path.
		for (int j = 1; j <= d; j++)
			average += path[j];
		average /= d;
		if (average > strike)
			return discount * (average - strike);
		else
			return 0.0;
	}

	public int getNumObsTimes() {
		return d;
	}


	public void simulate(RandomStream stream) {
		path = priceProcess.generatePath(stream);
	}


	public void simulateRuns(int n, RandomStream stream, Tally statValue,
			Tally statValuePos) {
		statValue.init();
		statValuePos.init();
		double x;
		for (int i = 0; i < n; i++) {
			simulate(stream);
			x = getPerformance();
			statValue.add(x);
			if (x > 0.0000000001)
				statValuePos.add(x);
			stream.resetNextSubstream();
		}
	}

	public String toString() {
		return "Asian option model with " + d + " observation times";
	}
	
	public AsianPriceProcess AsianPriceBGBS(double theta, double sigma, double nu, double s0, double r, double[] timeObs) {
		VarianceGammaProcess process_BGBS = new VarianceGammaProcess (0, 
				new BrownianMotion (0, theta, sigma, new MRG32k3a()), 
				new GammaProcessSymmetricalBridge(0, 1.0, nu, new MRG32k3a()));
		AsianPriceProcess asianPriceProcess_BGBS = new AsianPriceProcess(process_BGBS, s0, theta, sigma, r, timeObs);
		
		return asianPriceProcess_BGBS;
	}
	
	public AsianPriceProcess AsianPriceBGSS(double theta, double sigma, double nu, double s0, double r, double[] timeObs) {
		VarianceGammaProcess process_BGSS = new VarianceGammaProcess (0, 
				new BrownianMotion (0, theta, sigma, new MRG32k3a()), 
				new GammaProcess(0, 1.0, nu, new MRG32k3a()));
		AsianPriceProcess asianPriceProcess_BGSS = new AsianPriceProcess(process_BGSS, s0, theta, sigma, r, timeObs);
		
		return asianPriceProcess_BGSS;
	}
	
	public AsianPriceProcess AsianPriceDGBS(double theta, double sigma, double nu, double s0, double r, double[] timeObs) {
		VarianceGammaProcessDiff process_DGBS = new VarianceGammaProcessDiff(0, theta, sigma, nu, new MRG32k3a());
		AsianPriceProcess asianPriceProcess_DGBS = new AsianPriceProcess(process_DGBS, s0, theta, sigma, r, timeObs);
		
		return asianPriceProcess_DGBS;
	}
	
	
	public static void MCvsRQMCexperiment_SobolS(int sobolSize, int s, int n, int m, AsianOption asianOpt) {
		RandomStream stream = new MRG32k3a();
		
		DigitalNet p = new SobolSequence(sobolSize, 31, s); // n = 2^{16} points in s dim.
		PointSetRandomization rand = new RandomShift(new MRG32k3a());

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(asianOpt, stream, p, rand, n, m));
	}
	
	
	public static void MCvsRQMCexperiment_SobolLMS(int sobolSize, int s, int n, int m, AsianOption asianOpt) {
		RandomStream stream = new MRG32k3a();
		
		DigitalNet p = new SobolSequence(sobolSize, 31, s); // n = 2^{16} points in s dim.
		LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(asianOpt, stream, p, randLMS, n, m));
	}
	
	
	public static void MCvsRQMCexperiment_KorS(int sobolSize, int a, int s, int n, int m, AsianOption asianOpt) {
		RandomStream stream = new MRG32k3a();

		KorobovLattice pKor = new KorobovLattice (sobolSize, a, s);
		
		PointSetRandomization rand = new RandomShift(new MRG32k3a());

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(asianOpt, stream, pKor, rand, n, m));
	}
	
	
	public static void MCvsRQMCexperiment_KorBakerS(int sobolSize, int a, int s, int n, int m, AsianOption asianOpt) {
		RandomStream stream = new MRG32k3a();

		KorobovLattice pKor = new KorobovLattice (sobolSize, a, s);
		BakerTransformedPointSet pKorBaker = new BakerTransformedPointSet (pKor);
		
		PointSetRandomization rand = new RandomShift(new MRG32k3a());

		System.out.println (RQMCExperiment.makeComparisonExperimentMCvsRQMC 
				(asianOpt, stream, pKorBaker, rand, n, m));
	}

	
	public static void Derivate(RandomStream noise, int s, int d, double s0, double r, double theta, double sigma, double nu, 
								double delta, double[] timeObs, double K1) {
		
		GeometricVarianceGammaProcess asianPriceProcess_DGBS = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu, 
						   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise), 
						   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise)));
		
		GeometricVarianceGammaProcess asianPriceProcess_DGBSdelta = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu + delta, 
					   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise), 
					   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise)));
		
		AsianOption asianOpt_DGBS_d = new AsianOption(asianPriceProcess_DGBS, r, d, timeObs, K1);
		AsianOption asianOpt_DGBS_delta = new AsianOption(asianPriceProcess_DGBSdelta, r, d, timeObs, K1);

		int m = 100;
		int n = 1000;
		Tally statRQMC_sobol = new Tally ("RQMC averages for Asian option (with Sobol+LMS+S)");
		Tally statRQMC_Korobov = new Tally ("RQMC averages for Asian option (with Korobov+B+S)");
		Tally statDiff_CRN = new Tally("Stats on difference for MC with CRN");		
		Tally statDiff_IRN = new Tally("Stats on difference for MC with IRN");		
		statDiff_CRN.setConfidenceIntervalStudent();
		statDiff_IRN.setConfidenceIntervalStudent();
		statRQMC_sobol.setConfidenceIntervalStudent();
		statRQMC_Korobov.setConfidenceIntervalStudent();
		
		// RQMC WITH KOROBOV	
		
		int sobolSize = 16381;
		int a = 5693;
		KorobovLattice pKor = new KorobovLattice (sobolSize, a, s);
		BakerTransformedPointSet pKorBaker = new BakerTransformedPointSet (pKor);
		PointSetRandomization rand = new RandomShift(new MRG32k3a());
		RQMCExperiment.simulFDReplicatesRQMC (asianOpt_DGBS_d, asianOpt_DGBS_delta, delta, pKorBaker, rand, m,
				statRQMC_Korobov);
		System.out.println(statRQMC_Korobov.report(0.95, 6));
		
		
		// RQMC WITH SOBOL		
		DigitalNet pSobol = new SobolSequence(14, 31, s); 
		//LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());
		PointSetRandomization rand0 = new RandomShift(new MRG32k3a());
		RQMCExperiment.simulFDReplicatesRQMC (asianOpt_DGBS_d, asianOpt_DGBS_delta, delta, pSobol, rand0, m,
				statRQMC_sobol);
		
		System.out.println ("delta = " + delta);
		System.out.println(statRQMC_sobol.report(0.95, 6));
		System.out.println ("Variance per run: " + statRQMC_sobol.variance() * pSobol.getNumPoints() + "\n");
		
		// With simple MC CRN
		/*
		MonteCarloExperiment.simulFDReplicatesCRN (asianOpt_DGBS_d, asianOpt_DGBS_delta, delta, n, noise, statDiff_CRN);
	      System.out.println ("Ordinary MC with CRNs");
		  System.out.println(statDiff_CRN.report(0.95, 6));
	      System.out.println ("Variance per run: " + statDiff_CRN.variance() + "\n");
	      
	     MonteCarloExperiment.simulFDReplicatesIRN (asianOpt_DGBS_d, asianOpt_DGBS_delta, delta, n, noise, statDiff_IRN);
	      System.out.println ("Ordinary MC with IRNs");
		  System.out.println(statDiff_IRN.report(0.95, 6));
	      System.out.println ("Variance per run: " + statDiff_IRN.variance() + "\n");
	     */
	}
	
	
	public static class EquationQuestion_d implements MathFunction{
		
		double nu;
		double sigma ; 
		double K;
		double s0;
		double r;
		double theta;
		
		double omega;
		double mu_m, mu_p, nu_m, nu_p, lambda_m, lambda_p;
		
		double alpha;
		
		public EquationQuestion_d(double nu, double sigma, double K, double s0, double r, double theta) {
			this.nu = nu;
			this.sigma = sigma; 
			this.K = K;
			this.s0 = s0;
			this.r  = r;
			this.theta = theta;
			this.alpha = 1.0 / nu;
			
			this.omega = Math.log(1-(theta*nu)-(sigma*sigma*nu/2))/nu;
			
			this.mu_m = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) - theta);
			this.mu_p = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) + theta);
			this.nu_m = mu_m * mu_m * nu;
			this.nu_p = mu_p * mu_p * nu;
			
			this.lambda_m = mu_m / nu_m;
			this.lambda_p = mu_p / nu_p;
		}
		
		public double evaluate(double x){
			return  (alpha/(lambda_p - x) - alpha/(lambda_m + x)) - Math.log(K/s0) + r + omega;
		}
	}
	
	
	public static void main(String[] args) {
		
		double r  = 0.1;
		int d  = 8;
		double T  = 1;
		double nu = 0.2;
		double s0 = 100;
		double K1 = 90;
		double theta = -0.1436;
		double sigma = 0.12136;
		
		int rep = 32;   
		
		double[] timeObs = new double[d+1];
		for(int i=0; i<=d; i++) {
			timeObs[i] = (double)i/(double)d; 
			
		}
		
		
		//=========================================================================================================================================
		// PARTIE A
		//=========================================================================================================================================
		
		RandomStream noise = new MRG32k3a();
		
		// For RQMC experiment
		int s = 2*d;
		int sobolSize = 18;
		int a = 21876;
		int n = (int)Math.pow(2, 14);
		int m = 32;                     // Number of RQMC randomizations.

		
		
		System.out.println("=============================================================================================");
		System.out.println("BGSS");
		System.out.println("============================================================================================= \n");
		
		// BGSS: sequential sampling
		StochasticProcess asianPriceProcess_BGSS = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessAlternate(0.0,
				     new BrownianMotion (0.0, theta, sigma, new MRG32k3a()), 
				     new GammaProcess(0.0, 1.0, nu, noise)));

		AsianOption asianOpt_BGSS = new AsianOption(asianPriceProcess_BGSS, r, d, timeObs, K1);
		
		// Option value on 32 runs
		/*
		Tally statValuePos_BGSS = new Tally ("BGSS statValuePos");
		Tally statValue_BGSS = new Tally ("BGSS statValue");
		
		asianOpt_BGSS.simulateRuns(rep, new MRG32k3a(), statValue_BGSS, statValuePos_BGSS);
		System.out.println(statValue_BGSS.report());
		System.out.println(statValuePos_BGSS.report());
		*/
		
		
		//MCvsRQMCexperiment_SobolS(sobolSize, s, n, m, asianOpt_BGSS);
		//MCvsRQMCexperiment_SobolLMS(sobolSize, s, n, m, asianOpt_BGSS);
		//MCvsRQMCexperiment_KorS(sobolSize, a, s, n, m, asianOpt_BGSS);
		//MCvsRQMCexperiment_KorBakerS(sobolSize, a, s, n, m, asianOpt_BGSS);

		
		System.out.println("=============================================================================================");
		System.out.println("BGBS");
		System.out.println("============================================================================================= \n");
		
		// BGBS
		//GammaProcessSymmetricalBridge symetricGamma_BGBS = new GammaProcessSymmetricalBridge(s0, 1.0, nu, new MRG32k3a());
		//BrownianMotion bigBeautifulBM_BGBS = new BrownianMotion (s0, theta, sigma, new MRG32k3a());
		
		GeometricVarianceGammaProcess process_BGBS =  new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessAlternate(0.0,
						   new BrownianMotionBridge(0.0, theta, sigma, new MRG32k3a()), 
						   new GammaProcessSymmetricalBridge(0.0, 1.0, nu, noise)));

		AsianOption asianOpt_BGBS = new AsianOption(process_BGBS, r, d, timeObs, K1);
		
		// Option value on 32 runs
		/*
		Tally statValuePos = new Tally ("BGBS statValuePos");
		Tally statValue = new Tally ("BGBS statValue");
		
		asianOpt_BGBS.simulateRuns(rep, new MRG32k3a(), statValue, statValuePos);
		System.out.println(statValue.report());
		System.out.println(statValuePos.report());
		*/		
		
		//MCvsRQMCexperiment_SobolS(sobolSize, s, n, m, asianOpt_BGBS);
		//MCvsRQMCexperiment_SobolLMS(sobolSize, s, n, m, asianOpt_BGBS);
		//MCvsRQMCexperiment_KorS(sobolSize, a, s, n, m, asianOpt_BGBS);
		//MCvsRQMCexperiment_KorBakerS(sobolSize, a, s, n, m, asianOpt_BGBS);
		
		
		
		System.out.println("=============================================================================================");
		System.out.println("DGBS");
		System.out.println("============================================================================================= \n");
		
		// DGBS: difference sequential sampling
		GeometricVarianceGammaProcess asianPriceProcess_DGBS = new GeometricVarianceGammaProcess(s0, r,
				  new VarianceGammaProcessDiff(0.0,  theta,  sigma, nu, 
						   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise), 
						   new GammaProcessSymmetricalBridge(0.0, 1.0, 1.0, noise)));
		
		AsianOption asianOpt_DGBS = new AsianOption(asianPriceProcess_DGBS, r, d, timeObs, K1);
		
		// Option value on 32 runs
		/*
		Tally statValuePos_DBGS = new Tally ("DGBS statValuePos");
		Tally statValue_DBGS = new Tally ("DGBS statValue");
		
		asianOpt_DGBS.simulateRuns(rep, new MRG32k3a(), statValue_DBGS, statValuePos_DBGS);
		System.out.println(statValue_DBGS.report());
		System.out.println(statValuePos_DBGS.report());
		*/
		
		
		//MCvsRQMCexperiment_SobolS(sobolSize, s, n, m, asianOpt_DGBS);
		//MCvsRQMCexperiment_SobolLMS(sobolSize, s, n, m, asianOpt_DGBS);
		//MCvsRQMCexperiment_KorS(sobolSize, a, s, n, m, asianOpt_DGBS);
		//MCvsRQMCexperiment_KorBakerS(sobolSize, a, s, n, m, asianOpt_DGBS);
		

		
		//=========================================================================================================================================
		// PARTIE B
		//=========================================================================================================================================
		double delta = 0.001;
		Derivate(noise, s, d, s0, r, theta, sigma, nu, delta, timeObs, K1);
		
		
		
		//=========================================================================================================================================
		// PARTIE D
		//=========================================================================================================================================
		
		//double K = 130;
		//MathFunction equation = new EquationQuestion_d(nu, sigma, K, s0, r, theta);
		//double outil = RootFinder.brentDekker(0.0, 30, equation, 0.000000001);
		//System.out.println(outil);
	}
}
