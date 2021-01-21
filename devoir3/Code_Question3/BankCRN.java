package devoir3q3;
import java.io.IOException;

import umontreal.ssj.stat.Tally;

public class BankCRN extends Bank{
	
	   Tally statWait1_A = new Tally ("stats on derivate on the waiting time for A");
	   Tally statWait1_B = new Tally ("stats on derivate on the waiting for B");
	   
	   // 
	   double delta;

	   // The delta for the derivate
	   public BankCRN(double delta) throws IOException {
		   super();
		   this.delta = delta;
		}
	   
	   public void setS (double ss, double add) {
		      this.s = ss + add;  
	   }
	
	   // without crn
	   public void simulateDiff (int n) {
		  double w1_A, w1_B, w2_A, w2_B;
		  statWait1_A.init(); 
		  statWait1_B.init(); 
		  for (int i=0; i<n; i++) {
			 setS(10*MINUTE, 0);
		     simulateOneDay();
		     w1_A = statWaitsDay_A.sum() / nArrivals_A;
		     w1_B = statWaitsDay_B.sum() / nArrivals_B;
		     
		     setS(10*MINUTE, delta);
		     simulateOneDay();
		     w2_A = statWaitsDay_A.sum() / nArrivals_A;
		     w2_B = statWaitsDay_B.sum() / nArrivals_B;
		     statWait1_A.add ((w2_A - w1_A)/delta);
		     statWait1_B.add ((w2_B - w1_B)/delta);
		  }
	   }
	      
	  // with crn
	  public void simulateDiffCRN (int n) {
	      double w1_A, w1_B, w2_A, w2_B;
	      statWait1_A.init(); 
	      statWait1_B.init(); 
	
	      for (int i=0; i<n; i++) {
	         streamArr_A.resetNextSubstream();
	         streamServ_A.resetNextSubstream();
	         streamServ_B.resetNextSubstream();
	         streamIsClientBrdv.resetNextSubstream();
	         streamIsClientBcoming.resetNextSubstream();
	         streamDelayClientB.resetNextSubstream();
	         setS(10*MINUTE, 0);
	         simulateOneDay();
	         w1_A = statWaitsDay_A.sum() / nArrivals_A;
	         w1_B = statWaitsDay_B.sum() / nArrivals_B;
	
	
	         streamArr_A.resetNextSubstream();
	         streamServ_A.resetNextSubstream();
	         streamServ_B.resetNextSubstream();
	         streamIsClientBrdv.resetNextSubstream();
	         streamIsClientBcoming.resetNextSubstream();
	         streamDelayClientB.resetNextSubstream();
	         setS(10*MINUTE, delta);
	         simulateOneDay();
	         w2_A = statWaitsDay_A.sum() / nArrivals_A;
	         w2_B = statWaitsDay_B.sum() / nArrivals_B;
	         
	         statWait1_A.add ((w2_A - w1_A)/delta);
		     statWait1_B.add ((w2_B - w1_B)/delta);
	      }
	 }

	   public static void main (String[] args) throws IOException { 
		  BankCRN system = new BankCRN(10); 
		  
		  System.out.println ("For independant variable");
	      system.simulateDiff (5000);
	      System.out.println ("For client A");
	      system.statWait1_A.setConfidenceIntervalStudent();
	      System.out.println (system.statWait1_A.report (0.9, 3));
	      double varianceIndep_A = system.statWait1_A.variance();
	      System.out.println ("For client B");
	      system.statWait1_B.setConfidenceIntervalStudent();
	      System.out.println (system.statWait1_B.report (0.9, 3));
	      double varianceIndep_B = system.statWait1_B.variance();

	      System.out.println ("\n\n For CRN");
	      system.simulateDiffCRN (5000);
	      System.out.println ("For client A");
	      system.statWait1_A.setConfidenceIntervalStudent();
	      System.out.println (system.statWait1_A.report (0.9, 3));
	      double varianceCRN_A = system.statWait1_A.variance();
	      System.out.println ("For client B");
	      system.statWait1_B.setConfidenceIntervalStudent();
	      System.out.println (system.statWait1_B.report (0.9, 3));
	      double varianceCRN_B = system.statWait1_B.variance();

	      System.out.printf ("Variance ratio for client A:  %8.4g%n", varianceIndep_A/varianceCRN_A);
	      System.out.printf ("Variance ratio for client B:  %8.4g%n", varianceIndep_B/varianceCRN_B);
	      
	   }
}
