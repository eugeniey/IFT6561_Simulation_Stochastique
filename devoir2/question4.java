package devoir2;
import umontreal.ssj.probdist.LognormalDist;
import umontreal.ssj.rng.MRG32k3a;

public class question4 {
	public static void main(String[] args) {
		
		// Generator of uniform variables
		MRG32k3a uniformGenerator = new MRG32k3a();
		
		// Define constant
		double C = 100;
		double K = 102;
		double a = 100;
		double b = 102;
		double mu1 = 0.01;
		double mu2 = 0.01;
		double sigma1 = 0.04;
		double sigma2 = 0.04;
		
		double n = 10000;
		double confidence = 0.95;
		
		double sum_MC = 0.0;
		double sum_is = 0.0;
		double sum_MC_var = 0.0;
		double sum_is_var = 0.0;
		
		double X = 0.0;
		double X_is = 0.0;
		
		System.out.print("For b=102 \n"); 
		
		for(int i=0; i<=n; i++) {
			// Monte carlo
			double u1 = uniformGenerator.nextDouble();
			double w1 = LognormalDist.inverseF(mu1, sigma1, u1);
			double u2 = uniformGenerator.nextDouble();
			double w2 = LognormalDist.inverseF(mu2, sigma2, u2);
			
			double X_MC = 0.0;
			
			if(C*w1 <= a & C*w1*w2 >= b) {
				X_MC = C*w1*w2 - K;
			}
			sum_MC += X_MC;
			sum_MC_var += X_MC * X_MC;
			
			// Important sampling 
			double cdf1  = LognormalDist.cdf(mu1, sigma1, a/C);
			double w1_is = LognormalDist.inverseF(mu1, sigma1, u1 * cdf1);
			double cdf2  = LognormalDist.cdf(mu2, sigma2, b/(C*w1_is));
			double w2_is = LognormalDist.inverseF(mu2, sigma2, u2 * (1-cdf2) + cdf2);
			
			X = C*w1_is*w2_is - K;
			X_is = X * cdf1 * (1 - cdf2);
			
			sum_is += X_is;
			sum_is_var += X_is * X_is;
			
		}
		
		double average_MC = sum_MC/n;	
		double average_is = sum_is/n;
		double average_MC_square = sum_MC_var/n;	
		double average_is_square = sum_is_var/n;

		System.out.printf("Average Monte Carlo: %12.7f%n", average_MC);
		System.out.printf("Variance Monte Carlo: %12.7f%n", average_MC_square - (average_MC * average_MC));
		System.out.printf("Average important sampling: %12.7f%n", average_is);
		System.out.printf("Variance important sampling: %12.7f%n", average_is_square - (average_is * average_is));
		
		//========================================================================================================================
		
		double b2 = 112;
		
		sum_MC = 0;
		sum_MC_var = 0;
		sum_is = 0;
		sum_is_var = 0;
		
		for(int i=0; i<=n; i++) {
			// Monte carlo
			double u1 = uniformGenerator.nextDouble();
			double w1 = LognormalDist.inverseF(mu1, sigma1, u1);
			double u2 = uniformGenerator.nextDouble();
			double w2 = LognormalDist.inverseF(mu2, sigma2, u2);
			
			double X_MC = 0.0;
			
			if(C*w1 <= a & C*w1*w2 >= b2) {
				X_MC = C*w1*w2 - K;
			}
			sum_MC += X_MC;
			sum_MC_var += X_MC * X_MC;
			
			// Important sampling 
			double cdf1 = LognormalDist.cdf(mu1, sigma1, a/C);
			double w1_is = LognormalDist.inverseF(mu1, sigma1, u1 * cdf1);
			double cdf2 = LognormalDist.cdf(mu2, sigma2, b2/(C*w1_is));
			double w2_is = LognormalDist.inverseF(mu2, sigma2, u2 * (1-cdf2) + cdf2);
			
			X = C*w1_is*w2_is - K;
			X_is = X * cdf1 * (1 - cdf2);
			
			sum_is += X_is;
			sum_is_var += X_is * X_is;
		}
		
		average_MC = sum_MC/n;	
		average_is = sum_is/n;
		average_MC_square = sum_MC_var/n;	
		average_is_square = sum_is_var/n;
		
		System.out.print("\n\nFor b=112 \n"); 		
		System.out.printf("Average Monte Carlo: %12.7f%n", average_MC);
		System.out.printf("Variance Monte Carlo: %12.7f%n", average_MC_square - (average_MC*average_MC));
		System.out.printf("Average important sampling: %12.7f%n", average_is);
		System.out.printf("Variance important sampling: %12.7f%n", average_is_square - (average_is*average_is));				
	}
}

