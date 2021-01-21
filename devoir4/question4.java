package devoir4;

import java.util.LinkedList;
import java.util.ArrayList;
import java.util.Random;
import java.util.HashMap;
import java.util.Map;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import umontreal.ssj.probdist.ChiSquareDist;
import umontreal.ssj.rng.MRG32k3a;
import umontreal.ssj.rng.RandomStream; 

// Inspired from https://www.sanfoundry.com/java-program-graph-coloring-algorithm/


public class question4 {
	
	// number of colors
    private int nbColors;
    private String[] possibleColors = new String[]{"Red", "Green", "Blue", "Pink"};
    // number of nodes
    private int nbNode;
    // each index represents a node 
    // the value of color[i] is the color of the node i
    private String[] color; 
    // each index represents a node 
    // if the value of graph[i][j]=1 is mean nodes i and j are connected
    private int[][] graph;
    
    RandomStream stream2 = new MRG32k3a(); 
    
    
    public question4()
    {
    	nbColors = 4;
    	nbNode = 5;
    	graph = new int[nbNode][nbNode];
    	color = new String[nbNode];
    	
    	// fill the graph manually, not sure if it's the best idea
        for (int i = 0; i < nbNode; i++) {
        	for (int j = 0; j < nbNode; j++) {
            	graph[i][j] = 0;
            }
        }
        
        // node 1
    	graph[0][1] = 1;
    	graph[0][3] = 1;
    	graph[0][4] = 1;
    	// node 2
    	graph[1][0] = 1;
    	// node 3
    	graph[2][4] = 1;
    	// node 4
    	graph[3][0] = 1;
    	graph[3][4] = 1;
    	// node 5
    	graph[4][0] = 1;
    	graph[4][2] = 1;
    	graph[4][3] = 1;
    	
    	// Initiliase the color graph, how should we actually do this?
    	// Selon Gibbs: on part d'un état quelconque de l'ensemble d'état C
    	// we should be able to define it ourselves
    	color[0] = "Red"; //red
    	color[1] = "Pink"; //yellow
    	color[2] = "Pink"; //yellow
    	color[3] = "Blue"; //blue
    	color[4] = "Green"; //green
    }
    
    public String[] getColor() {
    	return this.color;
    }
    
    public String[] gibbs(String[] color) {
    	
    	// For each node
    	// we assign a new color randomly if possible
    	for(int currentnode=0; currentnode<nbNode; currentnode++) {
    		
    		LinkedList<String> illegalColor = new LinkedList<String>();
    		
    		// we look at all the neighbor
    		for(int neighbor=0; neighbor<nbNode; neighbor++) {
    			
        		// we look at his neighbor
    			// to look for the neighbors colors
        		if (graph[currentnode][neighbor] == 1) {
        			illegalColor.add(color[neighbor]);
        		}
        	}
    		
    		// After passing through the neighbors and looking at their colors
    		// we define the possible colors of the node
    		LinkedList<String> availableColor = new LinkedList<String>();
    		for(int c=0; c<nbColors; c++) {
    			if(!illegalColor.contains(possibleColors[c])) {
    				availableColor.add(possibleColors[c]);
    			}
        	}
    		
    		// We choose a random color in the available colors
    		int max = availableColor.size()-1;
    		int min = 0;
    		if(max > 0) {
    			Collections.shuffle(availableColor);
    			// this random is a bad choice, change it
    			Random rand = new Random();
    			int randomIndex = rand.nextInt(((max) - min) + 1) + min;
    			// assign the new color to the current node
    			color[currentnode] = availableColor.get(0);
    		}
    	}
    	
    	return color;
    }
    
    
    public static void main(String[] args) {
    	
    	// Map of the coloring and the counter of each coloring
    	Map<ArrayList<String>, Integer> coloringEncountered = new HashMap<ArrayList<String>, Integer>();
    	
    	question4 newGraph = new question4();
    	String[] colorState = newGraph.getColor();
    	
    	for(int i=0; i<216000; i++){

    		// We go into a new state
    		colorState = newGraph.gibbs(colorState);
    		
    		// cast as an array list for the key (cannot give an array as a key, bugs)
    		List<String> new_key_ = Arrays.asList(colorState);
    		ArrayList<String> new_key = new ArrayList<String>(new_key_);
    		
    		// If it's a coloring we already encountered
    		if(coloringEncountered.containsKey(new_key)) {
    			// update the counter
    			coloringEncountered.put(new_key, coloringEncountered.get(new_key) + 1);
    		}
    		// if it's the first time we get this coloring
    		else {    			
    			coloringEncountered.put(new_key, 1);
    		}
    	}
    	
    	double chi2_sum = 0;
    	double expected_obs = 216000/216;
    	
    	// print each results for fun
    	for (List<String> name: coloringEncountered.keySet()){
    		Integer num_observation = coloringEncountered.get(name);
            String value = num_observation.toString();  
            System.out.print("[");
            for(int i=0; i<name.size(); i++) {
            	System.out.print("node" + (i+1) + ":" + name.get(i) + " ");
            }
            System.out.print("]  ");  
            System.out.println(value);  
            
            chi2_sum += ((num_observation - expected_obs) * (num_observation - expected_obs))/ expected_obs;
            
    	} 
    	
    	RandomStream stream = new MRG32k3a(); 
    	// number of possibilities
    	System.out.println(coloringEncountered.keySet().size());  
    	System.out.println(chi2_sum);  
    	double fg = 1 - ChiSquareDist.cdf(215,25, chi2_sum);
    	System.out.println(fg); 
    	
    }
    
    
}
    		
