package travailFinal;

import travailFinal.AsianOption.EquationQuestion_d;
import umontreal.ssj.functions.MathFunction;
import umontreal.ssj.hups.DigitalNetBase2;
import umontreal.ssj.hups.LMScrambleShift;
import umontreal.ssj.hups.SobolSequence;
import umontreal.ssj.mcqmctools.MonteCarloExperiment;
import umontreal.ssj.mcqmctools.MonteCarloModelDouble;
import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;
import umontreal.ssj.stochprocess.BrownianMotion;
import umontreal.ssj.stochprocess.BrownianMotionBridge;
import umontreal.ssj.stochprocess.GammaProcess;
import umontreal.ssj.stochprocess.GammaProcessSymmetricalBridge;
import umontreal.ssj.stochprocess.GeometricVarianceGammaProcess;
import umontreal.ssj.stochprocess.StochasticProcess;
import umontreal.ssj.stochprocess.VarianceGammaProcessAlternate;
import umontreal.ssj.stochprocess.VarianceGammaProcessDiff;
import umontreal.ssj.util.RootFinder;


/*
 * Inspired from class exemples
 */
public class AsianOptionIS {
	
	double r  = 0.1;
	double T  = 1;
	double nu = 0.2;
	double s0 = 100;
	double theta = -0.1436;
	double sigma = 0.12136;
	
	double theta0 = 0.0;
	double K;
	
	int rep = 32;  
	
	int d = 1;   // One time step.
	
	double[] timeObs = new double[d+1];
	
	double discount;
   
    double dim = 2;
	double omega;

	double mu_m, mu_p; 
	double nu_m, nu_p;  
	double alpha;        
	double lambda_m, lambda_p;
	GammaDist gamma_m, gamma_p; 
	GammaDist gamma_m_is, gamma_p_is;  
	double L1;   
	double L2;      
    double condition;     
    
    Option option_regular = new Option ();
	OptionTwisting option_twisting = new OptionTwisting ();
	OptionTwistingCond option_twisting_cond = new OptionTwistingCond ();
	
	
	public AsianOptionIS (double K, double theta0) {
		
		this.K = K;
		this.theta0 = theta0;
		
		this.discount = Math.exp(-r * T);

		this.alpha = 1.0 / nu;
		this.omega = Math.log(1-(theta*nu)-(sigma*sigma*nu/2))/nu;

		this.mu_m = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) - theta);
		this.mu_p = 0.5 * (Math.sqrt(theta * theta + (2 * sigma * sigma / nu)) + theta);
		this.nu_m = mu_m * mu_m * nu;
		this.nu_p = mu_p * mu_p * nu;
		
		this.lambda_m = mu_m / nu_m;
		this.lambda_p = mu_p / nu_p;
		
		this.gamma_m = new GammaDist (alpha, lambda_m);
		this.gamma_p = new GammaDist (alpha, lambda_p);
		this.gamma_m_is = new GammaDist (alpha, lambda_m + theta0);
		this.gamma_p_is = new GammaDist (alpha, lambda_p - theta0);
		this.L1 = Math.pow((1.0 - theta0/lambda_p) * (1.0 + theta0/lambda_m), -alpha);
		this.L2 = Math.pow((1.0 + theta0/lambda_m), -alpha); 
		this.condition = Math.log(K/s0) - r - omega;
		
	}
	
	public class Option implements MonteCarloModelDouble {

		double value;

		public void simulate(RandomStream stream) {
			double gneg = gamma_m.inverseF(stream.nextDouble());
			double gpos = gamma_p.inverseF(stream.nextDouble());
			value = discount * (s0 * Math.exp(r + omega + gpos - gneg) - K);
			if (value < 0.0) value = 0.0;
		}

		public double getPerformance() {
			return value;
		}

    }
	

	public class OptionTwisting extends Option {

		public void simulate(RandomStream stream) {
			double gneg = gamma_m_is.inverseF(stream.nextDouble());
			double gpos = gamma_p_is.inverseF(stream.nextDouble());
			double L = Math.exp(theta0 * (gneg-gpos)) * L1;
			value = L * discount * (s0 * Math.exp(r + omega + gpos - gneg) - K); 
			if (value < 0.0) value = 0.0;			
		}
		public String toString() {
			return "European option on geometric VG process, one time step, IS twisting";
		}
	}
	
	
	public class OptionTwistingCond extends Option {

		public void simulate(RandomStream stream) {
			double gneg = gamma_m_is.inverseF(stream.nextDouble());
			double y = condition + gneg;   // We generate gpos conditional on > y
			double umin = gamma_p.cdf(y); 
			double u = umin + stream.nextDouble() * (1.0-umin);
			double gpos = gamma_p.inverseF(u);
			double L = L2 * Math.exp(theta0 * gneg) * 
					   gamma_p.barF(gneg + condition);
			value = (L * discount * (s0 * Math.exp(r + omega + gpos - gneg) - K));

			if (value < 0.0) value = 0.0;			
		}
	}
	

	
	public static class EquationQuestion_d implements MathFunction{
		
		double K;
		
		double r  = 0.1;
		double T  = 1;
		double nu = 0.2;
		double s0 = 100;
		double theta = -0.1436;
		double sigma = 0.12136;
		
		double omega;
		double mu_m, mu_p, nu_m, nu_p, lambda_m, lambda_p;
		
		double alpha;
		
		public EquationQuestion_d(double K) {
			this.K = K;
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
	
	public static double measureError(double var, double av) {
		double RMSE = Math.sqrt(var);
		return RMSE/av;
	}
	
	
	public static double optimise_theta3(double theta_init, double K, int n) {
		double min_p = 10000.0;
		double min_m = 10000.0;
		RandomStream noise = new MRG32k3a();
		Tally tally = new Tally("Estimate");
		
		double rate = 0.001;
		int search = 100000;
		
		double theta_p = theta_init;
		double theta_m = theta_init;
		
		double theta = theta_init;
		
		for(int i=0; i<search; i++) {
			
			if(theta>=37) {
				break;
			}
			
			AsianOptionIS test = new AsianOptionIS(K, theta);
			MonteCarloExperiment.simulateRunsDefaultReportStudent (test.option_twisting_cond, n, noise, tally, 0.95, 4);
			if(tally.variance() < min_p) {
				min_p = tally.variance();
				theta_p = theta;
			}
			
			theta += rate;
		}
		
		System.out.println("Done first");
		
		theta = theta_init;
		
		for(int i=0; i<search; i++) {
			
			if(theta<=0) {
				break;
			}
			
			AsianOptionIS test = new AsianOptionIS(K, theta);
			MonteCarloExperiment.simulateRunsDefaultReportStudent (test.option_twisting_cond, n, noise, tally, 0.95, 4);
			if(tally.variance() < min_m) {
				min_m = tally.variance();
				theta_m = theta;
			}
			
			theta -= rate;
		}
		
		if(min_m <= min_p) {
			return theta_m;
		}
		else {
			return theta_p;
		}	
	}
	
	
	public double getTheta_partE() {
		double Y = Math.log(this.K/this.s0) - this.r - this.omega;
		
		double sqrt = 4*this.alpha*this.alpha + this.lambda_m*this.lambda_m*Y*Y 
					  + 2*this.lambda_m*this.lambda_p*Y*Y + this.lambda_p*this.lambda_p*Y*Y;
		
		double numerateur = -2*this.alpha + Math.sqrt(sqrt) - this.lambda_m*Y + this.lambda_p*Y;

		return numerateur/(2*Y);
	}
	
	public static void main(String[] args) {
		
		RandomStream noise = new MRG32k3a();
		
		int dim = 2;
		
		
		//=========================================================================================================================================
		// PARTIE D
		//=========================================================================================================================================
		
		double K1 = 130;
		MathFunction equation1 = new EquationQuestion_d(K1);
		double outil1 = RootFinder.brentDekker(0.0, 30, equation1, 0.0000000000000001);
		//System.out.println(outil1);
		
		double K2 = 160;
		MathFunction equation2 = new EquationQuestion_d(K2);
		double outil2 = RootFinder.brentDekker(0.0, 26, equation2, 0.000000000000000001);

		
		int n = 100000;   // for MC.
		int m = 100;       // Number of replications.
		
		/*		
		AsianOptionIS test = new AsianOptionIS(K2, outil2);
		
		Tally statValue = new Tally("Stats on value of Asian option for MC");
		
		//System.out.println("  Ordinary MC:\n");
		//MonteCarloExperiment.simulateRunsDefaultReportStudent (test.option_regular, n, noise, statValue, 0.95, 4);
		//System.out.println(statValue.report());
		//System.out.println(measureError(statValue.variance(), statValue.average()));
		
		System.out.println("  MC with IS:\n");
		int nb_rea = 100;
		double[] avg_realisation = new double[nb_rea];
		
		for(int i=0; i<nb_rea; i++) {
			MonteCarloExperiment.simulateRunsDefaultReportStudent (test.option_twisting, n, noise, statValue, 0.95, 4);
			
			avg_realisation[i] = statValue.average();

		}
		
		double sum = 0;
		for (double d : avg_realisation) sum += d;

		double average_Xis = sum / nb_rea;
		
		double variance = 0;
		for (int i = 0; i < nb_rea; i++) {
		    variance += Math.pow(avg_realisation[i] - average_Xis, 2);
		}
		variance /= nb_rea;
		
		System.out.println("Average of average");
		System.out.println(average_Xis);
		System.out.print("Relative error: ");
		System.out.println(Math.sqrt(variance)/average_Xis);
		
		
		//System.out.print(statValue.report());
		//System.out.println(measureError(statValue.variance(), statValue.average()));
		

		// Then RQMC experiments.
		System.out.println("  RQMC:\n");
		DigitalNetBase2 pSobol = (new SobolSequence(16, 31, dim-1)).toNetShiftCj();  // 2^{16} points.
		LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());
		*/
		
		
		//=========================================================================================================================================
		// PARTIE E
		//=========================================================================================================================================
		AsianOptionIS test = new AsianOptionIS(K1, 1.0);
		double theta_from_e = test.getTheta_partE();
		
		AsianOptionIS test2 = new AsianOptionIS(K2, 1.0);
		double theta_from_e2 = test2.getTheta_partE();
		
		
		//=========================================================================================================================================
		// PARTIE F
		//=========================================================================================================================================
		
		AsianOptionIS asian_f = new AsianOptionIS(K2, theta_from_e2);
		
		Tally statValue = new Tally("Stats on value of Asian option for MC");
		
		System.out.println("  Ordinary MC:\n");
		MonteCarloExperiment.simulateRunsDefaultReportStudent (asian_f.option_twisting_cond, n, noise, statValue, 0.95, 4);
		System.out.println(statValue.report());

		/*
		int nb_rea = 100;
		double[] avg_realisation = new double[nb_rea];
		
		for(int i=0; i<nb_rea; i++) {
			MonteCarloExperiment.simulateRunsDefaultReportStudent (asian_f.option_twisting_cond, n, noise, statValue, 0.95, 4);
			
			avg_realisation[i] = statValue.average();

		}
		
		double sum = 0;
		for (double d : avg_realisation) sum += d;

		double average_Xis = sum / nb_rea;
		
		double variance = 0;
		for (int i = 0; i < nb_rea; i++) {
		    variance += Math.pow(avg_realisation[i] - average_Xis, 2);
		}
		variance /= nb_rea;
		
		System.out.println("Average of average");
		System.out.println(average_Xis);
		System.out.print("Relative error: ");
		System.out.println(Math.sqrt(variance)/average_Xis);
		*/
		
		//=================================
		// Optimisation
		//=================================
		
		/*
		Tally tally_IS_Cond = new Tally("Stats on value of Asian option for MC");
		
		//AsianOptionIS test = new AsianOptionIS(K1, 39);
	
		System.out.println("NUMERO f)");
		
		
		//System.out.println(outil2);
		//double opt = optimise_theta3(outil2, K2, 1000);
		//System.out.println(opt);
		AsianOptionIS asian_is_theta = new AsianOptionIS(K2, outil2);
		
		MonteCarloExperiment.simulateRunsDefaultReportStudent (asian_is_theta.option_twisting_cond, n, noise, tally_IS_Cond, 0.95, 20);
		
		System.out.print(tally_IS_Cond.report());
		System.out.print(tally_IS_Cond.variance());
		*/		
		
		//=========================================================================================================================================
		// PARTIE G
		//=========================================================================================================================================
		/*
		double theta_k1 = 22.574;
		double theta_k2 = 27.547;
		DigitalNetBase2 pSobol = (new SobolSequence(14, 31, 1)).toNetShiftCj();  // 2^{16} points.
		LMScrambleShift randLMS = new LMScrambleShift(new MRG32k3a());
		
		Tally tally_IS_Cond_sobol = new Tally("Stats on value of Asian option with IS COND for RQMC");
		
		tally_IS_Cond_sobol.setConfidenceIntervalStudent();
		
		AsianOptionIS asian_f = new AsianOptionIS(K2, theta_k2);
		TestOptionVGIS.simulReplicatesRQMCReport (asian_f.option_twisting_cond, pSobol, randLMS, m, tally_IS_Cond_sobol);
		*/
		
	}
}
