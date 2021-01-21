package devoir2;
import umontreal.ssj.rng.MRG32k3a;



public class question1 {
   public static void main(String[] args) {

      MathematicaSWB gen = new MathematicaSWB();
      MRG32k3a uniformGenerator = new MRG32k3a();
      
      double m = 10000;
      double ui = 0.0;
      double ui1 = 0.0;
      double ui2 = 0.0;
      
      // Nombre de collisions en moyenne pour les 3 options
      int C_1_avg = 0;
      int C_2_avg = 0;
      int C_3_avg = 0;
      
      for(int var = 0; var < 10; var++) {
    	  
    	  int C_1 = 0;
          int C_2 = 0;
          int C_3 = 0;
	      
          // définit le tableau de collisions pour les trois options
	      int[][][] collision_1 = new int[100][100][100];
	      int[][][] collision_2 = new int[100][100][100];
	      int[][][] collision_3 = new int[100][100][100];
	      double[] randomNumbers = new double[(int)(25*m+25)];
	      
	      // Initilisation collision array
	      for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				for (int k = 0; k < 100; k++) {
					collision_1[i][j][k] = 0;
				}
			}
	  	  }
	   // Initilisation collision array
	      for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				for (int k = 0; k < 100; k++) {
					collision_2[i][j][k] = 0;
				}
			}
	  	  }
	      // Initilisation collision array
	      for (int i = 0; i < 100; i++) {
			for (int j = 0; j < 100; j++) {
				for (int k = 0; k < 100; k++) {
					collision_3[i][j][k] = 0;
				}
			}
	  	  }
	      
	      // Initilisation random numbers 
	      for (int i = 0; i <= (int)(25*m+24); i++) {
	    	  	randomNumbers[i] = gen.nextDouble();
	  	  }
	      
	
	      for(int i=0; i<m; i++) {
		      ui  = Math.floor(randomNumbers[i]*100);
		      ui1 = Math.floor(randomNumbers[i+1]*100);
		      ui2 = Math.floor(randomNumbers[i+2]*100);
		      
		      collision_1[(int)ui][(int)ui1][(int)ui2] += 1;
		      
		      if(collision_1[(int)ui][(int)ui1][(int)ui2]>=2) {
					C_1 += 1;
			  }     
	      }
	      
	      for(int i=0; i<m; i++) {
		      ui  = Math.floor(randomNumbers[25*i]*100);
		      ui1 = Math.floor(randomNumbers[25*i + 20]*100);
		      ui2 = Math.floor(randomNumbers[25*i + 24]*100);
		      
		      collision_2[(int)ui][(int)ui1][(int)ui2] += 1;
		      
		      if(collision_2[(int)ui][(int)ui1][(int)ui2]>=2) {
					C_2 += 1;
			  }     
	      }
	      
	      for(int i=0; i<m; i++) {
		      ui  = Math.floor(uniformGenerator.nextDouble()*100);
		      ui1 = Math.floor(uniformGenerator.nextDouble()*100);
		      ui2 = Math.floor(uniformGenerator.nextDouble()*100);
		      
		      collision_3[(int)ui][(int)ui1][(int)ui2] += 1;
		      
		      if(collision_3[(int)ui][(int)ui1][(int)ui2]>=2) {
					C_3 += 1;
			  }     
	      }
	      
	      C_1_avg += C_1;
	      C_2_avg += C_2;
	      C_3_avg += C_3;
	      //System.out.print(C_3);
	      //System.out.print("\n");

      }
      
      System.out.print(C_1_avg/10);
      System.out.print("\n");
      System.out.print(C_2_avg/10);
      System.out.print("\n");
      System.out.print(C_3_avg/10);
      System.out.print("\n\n\n");
   }
}
