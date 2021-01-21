package devoir4;

import umontreal.ssj.probdist.GammaDist;
import umontreal.ssj.probdist.BetaDist;
import umontreal.ssj.probdist.NormalDist;
import umontreal.ssj.randvar.GammaGen;
import umontreal.ssj.randvar.NormalGen;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream;
import umontreal.ssj.stat.Tally;

public class question2 {
	
    double t0 = 0;
    double t1 = 0.6;
    double t2 = 0.8;
    double t3 = 1.0;
    double t4 = 1.2;
    
    double nu = 0.3;
    
    double alpha = 1/nu;
    double lambda = 1/nu;
    double theta = -0.1436;
    double mu = -0.1436;
    double sigma = 0.12136; 
	
	
	static Tally BGSS(int n, MRG32k3a uniformGen1, MRG32k3a uniformGen2) {
	    
		//double S = 0;
		double t0 = 0;
	    double t1 = 0.6;
	    double t2 = 0.8;
	    double t3 = 1.0;
	    double t4 = 1.2;
	    
	    double nu = 0.3;
	    double alpha = 1/nu;
	    double lambda = 1/nu;
	    double theta = -0.1436;
	    double mu = -0.1436;
	    double sigma = 0.12136; 
	    
	    double tau0 = 0;
	    double Yt0 = 0;
	    
	    Tally statBGSS = new Tally ("BGSS results");
	    statBGSS.init(); 
		
		for(int i=0; i<n; i++) {
			  double S = 0;
		      
		      // G(t1) - G(t0)
		      //double Gt1_Gt0 = GammaDist.inverseF ((t1-t0)*alpha, lambda, 15, uniformGen1.nextDouble());
			  double Gt1_Gt0 = GammaGen.nextDouble(uniformGen1, (t1-t0)*alpha, lambda);
		      double tau1 = Gt1_Gt0 + tau0;
		      // X(tau1) - X(tau0)
		      double Yt1 = NormalGen.nextDouble(uniformGen2, (tau1-tau0)*mu, (tau1-tau0)*Math.pow(sigma, 2)) + Yt0;
		      
		      // G(t2) - G(t1)
		      //double Gt2_Gt1 = GammaDist.inverseF ((t2-t1)*alpha, lambda, 15, uniformGen1.nextDouble());
		      double Gt2_Gt1 = GammaGen.nextDouble(uniformGen1, (t2-t1)*alpha, lambda);
		      double tau2 = Gt2_Gt1 + tau1;
		      // X(tau1) - X(tau0)
		      double Yt2 = NormalGen.nextDouble(uniformGen2, (tau2-tau1)*mu, (tau2-tau1)*Math.pow(sigma, 2)) + Yt1;
		      
		      // G(t3) - G(t2)
		      //double Gt3_Gt2 = GammaDist.inverseF ((t3-t2)*alpha, lambda, 15, uniformGen1.nextDouble());
		      double Gt3_Gt2 = GammaGen.nextDouble(uniformGen1, (t3-t2)*alpha, lambda);
		      double tau3 = Gt3_Gt2 + tau2;
		      // X(tau1) - X(tau0)
		      double Yt3 = NormalGen.nextDouble(uniformGen2, (tau3-tau2)*mu, (tau3-tau2)*Math.pow(sigma, 2)) + Yt2;
		      
		      // G(t4) - G(t3)
		      //double Gt4_Gt3 = GammaDist.inverseF ((t4-t3)*alpha, lambda, 15, uniformGen1.nextDouble());
		      double Gt4_Gt3 = GammaGen.nextDouble(uniformGen1, (t4-t3)*alpha, lambda);
		      double tau4 = Gt4_Gt3 + tau3;
		      // X(tau1) - X(tau0)
		      double Yt4 = NormalGen.nextDouble(uniformGen2, (tau4-tau3)*mu, (tau4-tau3)*Math.pow(sigma, 2)) + Yt3;

		      S = (Yt1 + Yt2 + Yt3 + Yt4);
		      
		      statBGSS.add(S);
	      }
   
	    return statBGSS;
	}
	
	
	static Tally BGBS(int n, MRG32k3a uniformGen1, MRG32k3a uniformGen2) {
		
		double t0 = 0;
	    double t1 = 0.6;
	    double t2 = 0.8;
	    double t3 = 1.0;
	    double t4 = 1.2;
	    
	    double nu = 0.3;
	    double alpha = 1/nu;
	    double lambda = 1/nu;
	    double theta = -0.1436;
	    double mu = -0.1436;
	    double sigma = 0.12136; 
	    
	    double Gt0 = 0;
	    double Yt0 = 0;
	    
	    Tally statBGBS = new Tally ("BGBS results");
	    statBGBS.init();
		
		for(int i=0; i<n; i++) {  
		  double S = 0;
			
		  double tau0 = 0;
		       
		  double tau4 = GammaDist.inverseF (t4*alpha, lambda, 15, uniformGen1.nextDouble());
		  double Gt4 = tau4;
		  //double tau4 = BetaDist.inverseF(G(ta), G(tb), uniformGen1.nextDouble());
		  // X(tau1) - X(tau0)
		  double Yt4 = NormalGen.nextDouble(uniformGen2, tau4*mu, tau4*Math.pow(sigma, 2));
		  
		  // beta
		  double beta2 = BetaDist.inverseF ((t2-t0)*alpha, (t4-t2)*alpha, uniformGen1.nextDouble());
		  double tau2 = beta2 * (Gt4 - Gt0) + Gt0;
		  // X(tau1) - X(tau0)
		  double Yt2 = NormalGen.nextDouble(uniformGen2, Yt0 + (Yt4 - Yt0)*(tau2/tau4), Math.pow(sigma, 2)*tau2*(tau4-tau2)/tau4);
		  
		  // G(t1) - G(t0)
		  double beta1 = BetaDist.inverseF ((t1-t0)*alpha, (t2-t1)*alpha, uniformGen1.nextDouble());
		  double tau1 = beta1 * (tau2 - tau0) + tau0;
		  // X(tau1) - X(tau0)
		  double Yt1 = NormalGen.nextDouble(uniformGen2, Yt0 + (Yt2 - Yt0)*(tau1/tau2), Math.pow(sigma, 2)*tau1*(tau2-tau1)/tau2);
		
		  // G(t3) - G(t2)
		  double beta3 = BetaDist.inverseF ((t3-t1)*alpha, (t4-t3)*alpha, uniformGen1.nextDouble());
		  double tau3 = beta3 * (tau4 - tau2) + tau2;
		  // X(tau1) - X(tau0)
		  double Yt3 = NormalGen.nextDouble(uniformGen2, Yt2 + (Yt4 - Yt2)*(tau3/tau4), Math.pow(sigma, 2)*tau3*(tau4-tau3)/tau4);
		
		  
		  S = (Yt1 + Yt2 + Yt3 + Yt4);
		  statBGBS.add(S); 
		}

		return statBGBS;
	}
	
	
	static Tally DGBS(int n, MRG32k3a uniformGen1, MRG32k3a uniformGen2) {
		
		double t0 = 0;
	    double t1 = 0.6;
	    double t2 = 0.8;
	    double t3 = 1.0;
	    double t4 = 1.2;
	    
	    double nu = 0.3;
	    double alpha = 1/nu;
	    double lambda = 1/nu;
	    double theta = -0.1436;
	    double mu = -0.1436;
	    double sigma = 0.12136; 
	    
	    double Gt0 = 0;
	    double Yt0 = 0;
	    
	    double mu_p = (Math.sqrt(theta*theta + 2*sigma*sigma/nu) + theta) / 2;
	    double mu_m = (Math.sqrt(theta*theta + 2*sigma*sigma/nu) - theta) / 2;
	    double nu_p = mu_p*mu_p*nu; 
	    double nu_m = mu_m*mu_m*nu;
	    
	    Tally statDGBS = new Tally ("DGBS results");
	    statDGBS.init();
		
		for(int i=0; i<n; i++) {  
		  double S = 0;
			
		  double tau0 = 0;
		  
		  double alpha_p = mu_p*mu_p/nu_p;
		  double alpha_m = mu_m*mu_m/nu_m;
		       
		  double G4_p = GammaDist.inverseF (t4*alpha_p, mu_p/nu_p, 15, uniformGen1.nextDouble());
		  double G4_m = GammaDist.inverseF (t4*alpha_m, mu_m/nu_m, 15, uniformGen2.nextDouble());
		  //double tau4 = BetaDist.inverseF(G(ta), G(tb), uniformGen1.nextDouble());
		  // X(tau1) - X(tau0)
		  double Yt4 = G4_p - G4_m;
		  
		  // beta
		  double beta2_p = BetaDist.inverseF ((t2-t0)*alpha_p, (t4-t2)*alpha_p, uniformGen1.nextDouble());
		  double beta2_m = BetaDist.inverseF ((t2-t0)*alpha_m, (t4-t2)*alpha_m, uniformGen2.nextDouble());
		  double G2_p = beta2_p * (G4_p - Gt0) + Gt0;
		  double G2_m = beta2_m * (G4_m - Gt0) + Gt0;
		  // X(tau1) - X(tau0)
		  double Yt2 = G2_p - G2_m;
		  
		  // G(t1) - G(t0)
		  double beta1_p = BetaDist.inverseF ((t1-t0)*alpha_p, (t2-t1)*alpha_p, uniformGen1.nextDouble());
		  double beta1_m = BetaDist.inverseF ((t1-t0)*alpha_m, (t2-t1)*alpha_m, uniformGen2.nextDouble());
		  double G1_p = beta1_p * (G2_p - Gt0) + Gt0;
		  double G1_m = beta1_m * (G2_m - Gt0) + Gt0;
		  // X(tau1) - X(tau0)
		  double Yt1 = G1_p - G1_m;
		
		  // G(t3) - G(t2)
		  double beta3_p = BetaDist.inverseF ((t3-t2)*alpha_p, (t4-t3)*alpha_p, uniformGen1.nextDouble());
		  double beta3_m = BetaDist.inverseF ((t3-t2)*alpha_m, (t4-t3)*alpha_p, uniformGen2.nextDouble());
		  double G3_p = beta3_p * (G4_p - G2_p) + G2_p;
		  double G3_m = beta3_m * (G4_m - G2_m) + G2_m;
		  // X(tau1) - X(tau0)
		  double Yt3 = G3_p - G3_m;
		
		  
		  S = (Yt1 + Yt2 + Yt3 + Yt4);
		  statDGBS.add(S); 
		}

		return statDGBS;
	}
	
	
	public static void main(String[] args) {

	      MRG32k3a uniformGen1  = new MRG32k3a();
	      MRG32k3a uniformGen2 = new MRG32k3a();
	      MRG32k3a uniformGen1BGBS  = new MRG32k3a();
	      MRG32k3a uniformGen2BGBS = new MRG32k3a();
	      MRG32k3a uniformGen1DGBS  = new MRG32k3a();
	      MRG32k3a uniformGen2DGBS = new MRG32k3a();

	      int n = 10000;
	      
	      Tally bgssTally = BGSS(n, uniformGen1, uniformGen2);
	      Tally bgbsTally = BGBS(n, uniformGen1BGBS, uniformGen2BGBS);
	      Tally dgbsTally = DGBS(n, uniformGen1DGBS, uniformGen2DGBS);
	      
	      System.out.println(bgssTally.report());
	      System.out.println(bgssTally.formatCINormal(0.95, 4));

	      System.out.println(bgbsTally.report());
	      System.out.println(bgbsTally.formatCINormal(0.95, 4));
	      
	      System.out.println(dgbsTally.report());
	      System.out.println(dgbsTally.formatCINormal(0.95, 4));
	   }
}
