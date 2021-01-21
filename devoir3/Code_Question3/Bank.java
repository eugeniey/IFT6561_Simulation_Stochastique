package devoir3q3;

import umontreal.ssj.simevents.*;
import umontreal.ssj.rng.*;
import umontreal.ssj.randvar.*;
import umontreal.ssj.probdist.*;
import umontreal.ssj.stat.Tally;
import java.io.*;
import java.util.StringTokenizer;


//import devoir3.CallCenter.Call;
//import devoir3.CallCenter.CallCompletion;
//import devoir3.CallCenter.NextPeriod;

import java.util.LinkedList;

public class Bank {
	
   static final double HOUR = 3600.0;  // Time is in seconds. 
   static final int MINUTE = 60;  // Time is in seconds. 

   // Data
   // Arrival rates are per hour, service and patience times are in seconds.
   int openingTime;    // Opening time of the center (in hours).
   int numPeriods;     // Number of working periods (hours) in the day.
   
   int[] numAgents_C;  // Number of agents for each period. 
   int[] numAgents_D;  // Number of agents for each period. 
   
   
   double r;              // probability that the plage has a rdv
   
   // Info for client A
   double[] lambda;     // Base arrival rate lambda_j for each j for client A 
   double mu_A;           // Average of the lognormal for client A 
   double sigma_A;        // Écart type of the lognormal for client A 
   // Wait list and arrive for client A
   LinkedList<Client_A> waitList_A  = new LinkedList<Client_A>();
   Event nextArrival = new Arrival();   

   
   // Info for client B
   double alpha0;         // Parameter of gamma distribution for B.
   double p;              // Probability that patience time is 0.
   double nu;             // Parameter of exponential for patience time.
   double alpha, gamma;   // Parameters of gamma service time distribution. 
   double s;              // Want stats on waiting times smaller than s.
   int currentPeriod;
   double mu_B;           // Average of the service time  for client B
   double sigma_B;        // Écart type of the service time for client B
   double mu_B_r;         // Average to define lateness for client B
   double sigma_B_r;      // Écart type to define lateness client B
   LinkedList<LinkedList<Client_B>>[] clientsOfAgents_DforTheDay; // list of agent D and their clients B for the day
   Service_B processForClient_B;  // Processus de rendez-vous pour client B
   
   // Variables
   double arrRate = 0.0;    // Current arrival rate.
   int nAgents_C;           // Number of agents C in current period.
   int nAgents_D;           // Number of agents D in current period.
   int nBusy_C;             // Number of agents C occupied;
   int nArrivals_A;         // Number of arrivals today for client A;
   int nArrivals_B;         // Number of arrivals today for client B;

   RandomStream streamArr_A = new MRG32k3a();           // For arrivals A
   RandomStream streamServ_A = new MRG32k3a();          // For service A
   RandomStream streamServ_B = new MRG32k3a();          // For service B
   RandomStream streamIsClientBrdv = new MRG32k3a();    // If client B has rdv
   RandomStream streamIsClientBcoming = new MRG32k3a(); // If client B shows up
   RandomStream streamDelayClientB = new MRG32k3a();    // Define lateness of client B
   
   // Define Tallies for different statistics
   Tally statArrivals_A = new Tally ("Number of arrivals per day");
   Tally statWaitsDay_A = new Tally ("Waiting times within a day for client A");
   Tally statWaitsDay_B = new Tally ("Waiting times within a day for client B");
   Tally statWaits_A    = new Tally ("Average waiting time per client A");
   Tally statWaits_B    = new Tally ("Average waiting time per client B");
 
   
   public Bank() throws IOException {
	   
	  // Define the parameters
      alpha0     = 10.0;
      p          =  0.05 ;
      r          = 0.8;
      nu         =  0.001 ;
      alpha      = 1.0 ;  
      gamma      = 0.01 ;
      s          = 10 * MINUTE;
      numPeriods = 3;
      mu_A       = 240;
      sigma_A    = 50;
      mu_B       = 25 * MINUTE;
      sigma_B    = 4 * MINUTE;
      mu_B_r     = -100;
      sigma_B_r  = 60;
      
      numAgents_C = new int[] {3,4,3}; 
      numAgents_D = new int[] {2,3,2}; 
      lambda = new double[] {30,40,25};         
   }
   

   // Event: A client type A arrives.
   class Arrival extends Event {
      public void actions() {
         nextArrival.schedule 
            (ExponentialDist.inverseF (arrRate, streamArr_A.nextDouble()));
         nArrivals_A++;
         new Client_A();               // Call just arrived.
      }
   }
   
   
   /*
    * Define information on client A object
    */
   class Client_A {
	   double arrivalTime; 
	   double serviceTime;
	   
	   public Client_A() {
		   // Generate service time.
		   LognormalDistFromMoments logNormal = new LognormalDistFromMoments(mu_A, Math.pow(sigma_A, 2));
		   this.serviceTime = logNormal.inverseF(streamServ_A.nextDouble());
		   // Start service immediately.
		   if (nBusy_C < nAgents_C) {           
	        	 nBusy_C++;
	            statWaitsDay_A.add(0.0);
	            new ClientCompletion_A().schedule(serviceTime);
	         } 
	         else {
	        	 arrivalTime = Sim.time();
		         waitList_A.addLast(this); 
	         }
	   }
	   
	   public void endWait(boolean isAgentD) {
		     if(!isAgentD) nBusy_C++;
	         double wait = Sim.time() - arrivalTime;
	         new ClientCompletion_A().schedule(serviceTime);
	         statWaitsDay_A.add (wait);
	    }
   }
     
   // Event: A call is completed.
   class ClientCompletion_A extends Event {
      public void actions() { nBusy_C--;   checkQueue(); }
   }  
   
   // Start answering new calls if agents are free and queue not empty.
   public void checkQueue() {
      while ((waitList_A.size() > 0) && (nBusy_C < nAgents_C))
         ((Client_A)waitList_A.removeFirst()).endWait(false);
   }
   
     
   /*
    * Define information on client B object
    */
   class Client_B {
	   double arrivalTime;    // time client arrive
	   double serviceTime;    // time of service
	   double endServiceTime; // time at which the service ends
	   boolean isComing;      // if the clients shows up
	   boolean isRDV;         // if the client has a rdv
	   double delay;          // the delay of the client (can be negative)
	   double appoinTime;     // time of the appointment
	   double startServiceTime;  // time at which the service started  
	   
	   // Define information on Client B
	   public Client_B(boolean isComing, boolean isRDV, double appoinTime, double delay) {
		   this.isComing    = isComing;
		   this.appoinTime  = appoinTime;
		   this.delay       = delay;
		   this.isRDV       = isRDV;
		   if(delay<=0) {
			   this.arrivalTime = appoinTime;
		   }
		   else {
			   this.arrivalTime = appoinTime + delay;
		   }
		   
		   
		   LognormalDistFromMoments logNormal = new LognormalDistFromMoments(mu_B, Math.pow(sigma_B, 2));
		   this.serviceTime = logNormal.inverseF(streamServ_B.nextDouble());
	   }

	   // Define the end of service and wait time for the client
	   public void defineEndServiceTimeAndWait() {
		   this.endServiceTime = Sim.time();
		   // We just add a waiting time is the client shows up or has a rdv
		   if(this.isComing && this.isRDV) {
			   double startServiceTime = Sim.time() - this.serviceTime;
			   double wait = startServiceTime - arrivalTime;
			   // if client comes before being served
			   if(wait>0) {
				   statWaitsDay_B.add(wait); 
			   }
			   else {
				   statWaitsDay_B.add(0.0); 
			   }
		   }
	   }
   }
   

   /*
    * Generate the service  for the whole day 
    * for each agent D and each of their clients B
    */
   class Service_B extends Event {

	   Client_B currentClient;      // the current client
	   int indexCurrentClient;      // define the current client
	   int indexAgent;              // define the agent
	   int period;					// the period at which it was instanciated
	   int nbPeriodAgentIsThere;    // define for how many period the agent is there
	   boolean clientThatPassAgain; // if client didn't get his service done
	   
	   // Define information for service B.    
	   public Service_B(int indexAgent, int period) {
		   this.period = period;
		   this.indexAgent = indexAgent;
		   // the current client is always the first one of the list of the agent
		   this.indexCurrentClient = 0;
		   if(period == 0) {
			   this.nbPeriodAgentIsThere = 2;
		   }
		   else {
			   this.nbPeriodAgentIsThere = 1;
		   }
	   }

	   // Schedule at which point the service will be done for the current client
	   // and will call the next client
	   public void nextClient() {

		   // Here Sim.time() is the end of the previous rdv
		   if(indexCurrentClient <= 3) {   
			   
			   // get the current client 
			   currentClient = clientsOfAgents_DforTheDay[this.period].get(indexAgent).get(indexCurrentClient);
			   
			   // End of appointment in case a client doesn't show up
			   double endOfService = currentClient.appoinTime + 30*MINUTE;
			   
			   // If there is an appointment
			   if(currentClient.isRDV) {
				   // Time before the appointment
				   double timeBeforeNextAppoint = currentClient.appoinTime - Sim.time();
				   // Time before the client actually comes
				   double timeBeforeNextArrival = currentClient.arrivalTime - Sim.time();

				   // If there is more than s seconds before the next appointment
				   if(timeBeforeNextAppoint >= s) {
					   // If there is waiting client A
					   if(waitList_A.size() > 0) {
						   Client_A servedClientA = (Client_A)waitList_A.removeFirst();
						   servedClientA.endWait(true);  // Serve a type A client
						   
						   double timeDoneWithClientA = Sim.time() + servedClientA.serviceTime;
						   
						   // After he finished with client A, agent D goes to his appointment
						   // If client B comes before agent D is done with client A
						   if(timeDoneWithClientA > currentClient.arrivalTime) {
							   // He answers his next client when clientA is done
							   if(currentClient.isComing) {
								   this.schedule(servedClientA.serviceTime + currentClient.serviceTime);
							   }
							   // If client doesn't come, we schedule end of service at the next meeting
							   // measure the time left until end of appointment
							   else {
								   if(endOfService - timeDoneWithClientA < 0) {
									   this.schedule(servedClientA.serviceTime);  
								   }
								   else {
									   this.schedule(endOfService - timeDoneWithClientA);   
								   }
							   }
						   }
						   // If client B comes after agentD is done with client A
						   else {
							   // He answers his next client when clientA is done
							   if(currentClient.isComing) {
								   // timeBeforeNextArrival cannot be smaller than 0
								   this.schedule(Math.abs(timeBeforeNextArrival) + currentClient.serviceTime);
							   }
							   // If client doesn't come, we schedule end of service at the next meeting
							   // measure the time left until end of appointment
							   else {
								   this.schedule(endOfService - timeDoneWithClientA);
							   }
						   } 
						}
					   
					   // No client A to serve, he waits until clients B comes
					   else {
						   if(currentClient.isComing) {
							   // If person comes before previous appoint is done
							   if(timeBeforeNextArrival <= 0) {
								   this.schedule(currentClient.serviceTime);
							   }
							   // If person comes after previous appoint is done
							   else {
								   this.schedule(Math.abs(timeBeforeNextArrival) + currentClient.serviceTime);
							   }
						   }
						   // If client doesn't come, we schedule end of service at the next meeting
						   else {
							   this.schedule(endOfService - Sim.time());
						   }
					   }
				   }
				   // If there is less than s second before appointment
				   // agentD has to wait until client B comes
				   // End of current service is scheduled based on the arrival and previous end schedule
				   else{
					   if(currentClient.isComing) {
						   // If person comes before previous appoint is done
						   if(timeBeforeNextArrival <= 0) {
							   this.schedule(currentClient.serviceTime);
						   }
						   // If person comes after previous appoint is done
						   else {
							   this.schedule(Math.abs(timeBeforeNextArrival) + currentClient.serviceTime);
						   }
					   }
					   // If client doesn't come, we schedule end of service at the next meeting
					   else {
						   this.schedule(endOfService - Sim.time());
					   }
				   }				   
			   }
			   // If agent D doesn't have an appointment
			   // can complete client A service
			   else {
				   // If there is waiting client A
				   if(waitList_A.size() > 0) {
					   Client_A servedClientA = (Client_A)waitList_A.removeFirst();
					   servedClientA.endWait(true);  // Serve a type A client
					   
					   // After he finished with client A, agent D had nothing to do
					   // Consider end time when done with client A
					   this.schedule(servedClientA.serviceTime);
					}
				   // If there is no client A to serve, waits until next client
				   else {
					   this.schedule(endOfService - Sim.time());
				   }
			   }
			   
			   // prepare for next client
			   indexCurrentClient++; 
		   }
		   
		   // We go into a new period
		   else if(this.period < nbPeriodAgentIsThere) {
			   this.period = this.period + 1;
			   indexCurrentClient = 0;
			   // the client has to come again
			   this.clientThatPassAgain = true;
			   
			   // this means that the current client is one from a new period
			   // we need to make it pass into the nextClient again, but with the new information
			   this.schedule(0.0);
		   }
	   }
	   
	   public void actions() {
		   // Define end service time only if it's not a empty client 
		   if(!this.clientThatPassAgain) {
			   currentClient.defineEndServiceTimeAndWait();	 
		   }
		   else {
			   this.clientThatPassAgain = false;
		   }
		   this.nextClient(); 
	   }
   }
   

   
   // Event: A new period begins.
   class NextPeriod extends Event {
      int j;     // Number of the new period.

      public NextPeriod(int period){ 
    	  j = period; 
    	  currentPeriod = j;
      }
      public void actions() {

         if (j < numPeriods) {
        	nAgents_C = numAgents_C[j];
            arrRate = (double)lambda[j]/(double)HOUR;
           
            if (j == 0) {
            	// We start the service of each agent 
                for(int a=0; a<numAgents_D[j]; a++) {
                	processForClient_B = new Service_B(a, j);
                	processForClient_B.nextClient();
                }
            	nextArrival.schedule (ExponentialDist.inverseF (arrRate, streamArr_A.nextDouble()));
            }
               
            else {
            	// if we get an agent more in a certain period
            	if(numAgents_D[j] > numAgents_D[j-1]) {
	        		for(int a=numAgents_D[j-1]; a<numAgents_D[j]; a++) {
	                	processForClient_B = new Service_B(a, j);
	                	processForClient_B.nextClient();
	                }
            	}
               checkQueue();
               nextArrival.reschedule ((nextArrival.time() - Sim.time()) * lambda[j-1] / lambda[j]);
            }
            new NextPeriod(j+1).schedule (2.0 * HOUR);
         }
         else
            nextArrival.cancel();  // End of the day.
      }
   }
   
   
   public void simulateOneDay() { 
	      Sim.init();   
	      statWaitsDay_A.init();
	      statWaitsDay_B.init();
	      nArrivals_A = 0;   
	      nArrivals_B = 0;    
	      nBusy_C = 0;
	      
	      clientsOfAgents_DforTheDay = new LinkedList[numPeriods];
	      
	      // Initialise all the agents and their clients for each period
	      for(int i=0; i<numPeriods; i++) {
	    	  
	    	  // List of agents D for the period
	    	  LinkedList<LinkedList<Client_B>> allClientsOfAgents_D = new LinkedList<LinkedList<Client_B>>();
	    	  
	    	  for(int k=0; k<numAgents_D[i]; k++) {

	    		  // List of clients B for the agent
		    	  LinkedList<Client_B> clients_BForOnePeriod = new LinkedList<Client_B>();
		    	  
		    	  // Create the clients B for the period
		    	  for(int j=0; j<4; j++) {
		    		  
		    		  double plannedAppoinTime = openingTime + j*30*MINUTE + i*2*HOUR;
		    		  
		    		  // If there is a rdv 
		    		  if(streamIsClientBrdv.nextDouble() < r) {
		    			  
			    		  // Client doesn't show up
			    		  if(streamIsClientBcoming.nextDouble() < p) {
			    			  Client_B client_notcoming = new Client_B(false, true, plannedAppoinTime, 0); 
				    		  clients_BForOnePeriod.addLast(client_notcoming);
			    		  }
			    		  // Client shows up
			    		  else {
			    			  Client_B client_coming;
			    			  if(j == 0) {
			    				  double delay = NormalDist.inverseF(mu_B_r, sigma_B_r, streamDelayClientB.nextDouble());
			    				  // Because client cannot be there before opening. 
			    				  if(delay>0) {
			    					  client_coming = new Client_B(true, true, plannedAppoinTime, delay); 
			    				  }else {
			    					  client_coming = new Client_B(true, true, plannedAppoinTime, 0.0); 
			    				  }
			    			  }
			    			  // Clients can be early
			    			  else {
			    				  double delay = NormalDist.inverseF(mu_B_r, sigma_B_r, streamDelayClientB.nextDouble());
				    			  client_coming = new Client_B(true, true, plannedAppoinTime, delay); 
			    			  }
			    			  
				    		  clients_BForOnePeriod.addLast(client_coming);
				    		  nArrivals_B++;
			    		  }
		    		  }
		    		  // If there is not rdv planned
		    		  else {
		    			  Client_B client_noRDV = new Client_B(false, false, plannedAppoinTime, 0); 
		    			  clients_BForOnePeriod.addLast(client_noRDV);
		    		  } 
		    	  }
		    	  allClientsOfAgents_D.addLast(clients_BForOnePeriod);
		    	  
	    	  }
	    	  clientsOfAgents_DforTheDay[i] = allClientsOfAgents_D;
	      }
	      
	      
	      new NextPeriod(0).schedule(openingTime * HOUR);
	      Sim.start();
	      // Here the simulation is running...
	      statArrivals_A.add ((double)nArrivals_A);      
	      statWaits_A.add (statWaitsDay_A.sum() / nArrivals_A);
	      statWaits_B.add (statWaitsDay_B.sum() / nArrivals_B);
	      
	      System.out.println(statWaitsDay_B.sum() / nArrivals_B);
   }

   
   static public void main (String[] args) throws IOException { 

	   	  // Instanciation a bank
		  Bank bank = new Bank() ; 
		  
		  // simulate 1000 days
	      for (int i=1; i <= 1000; i++) {
	    	  bank.simulateOneDay();
	      }
	      
	      System.out.println (
	        bank.statArrivals_A.reportAndCIStudent (0.9) +
	        "\n" +
	        bank.statWaits_A.reportAndCIStudent (0.9) +
	        "\n" +
	        bank.statWaits_B.reportAndCIStudent (0.9) ); 
   }


}
