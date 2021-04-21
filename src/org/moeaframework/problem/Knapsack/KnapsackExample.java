package org.moeaframework.problem.Knapsack;

/* Copyright 2009-2016 David Hadka
*
* This file is part of the MOEA Framework.
*
* The MOEA Framework is free software: you can redistribute it and/or modify
* it under the terms of the GNU Lesser General Public License as published by
* the Free Software Foundation, either version 3 of the License, or (at your
* option) any later version.
*
* The MOEA Framework is distributed in the hope that it will be useful, but
* WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
* or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
* License for more details.
*
* You should have received a copy of the GNU Lesser General Public License
* along with the MOEA Framework.  If not, see <http://www.gnu.org/licenses/>.
*/


import java.io.IOException;
import java.io.InputStream;
import java.util.ArrayList;
import java.util.Vector;
import java.util.regex.Matcher;
import java.util.regex.Pattern;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.Reader;
import org.moeaframework.Executor;
import org.moeaframework.util.*;
import org.moeaframework.analysis.plot.Plot;
import org.moeaframework.Analyzer;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.Knapsack.Knapsack;
import org.moeaframework.util.ReferenceSetMerger;
import org.moeaframework.util.io.CommentedLineReader;
import org.moeaframework.core.PopulationIO;

/**
 * Example of binary optimization using the {@link Knapsack} problem on the
 * {@code knapsack.100.2} instance.
 */
public class KnapsackExample {
	
	private int nsacks;
	private int nitems;
	private int nobjectives;
	private int[][] weight;
	private int[][] profit;
	private int[] capacity;
	private String[] itemNames;
	
	public KnapsackExample(String path)throws IOException {
		 load (path);
	}
	
	private void load(String path) throws IOException {
		 Reader reader = new FileReader(path);
		 Pattern specificationLine = Pattern.compile("knapsack problem specification \\((\\d+) knapsacks, (\\d+) items, (\\d+) objectives\\)");
		 Pattern capacityLine = Pattern.compile(" capacity: \\+(\\d+)");
		 Pattern weightLine = Pattern.compile("  weight: \\+(\\d+)");
		 Pattern profitLine = Pattern.compile("  profit: \\+(\\d+)");
		 Pattern itemLine = Pattern.compile("//item (\\d+):(...)");

		CommentedLineReader lineReader = null;
		String line = null;
		Matcher matcher = null;
		
		try {
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); // the problem specification line
			matcher = specificationLine.matcher(line);
				
			if (matcher.matches()) {
				nsacks = Integer.parseInt(matcher.group(1));
				nitems = Integer.parseInt(matcher.group(2));
				nobjectives = Integer.parseInt(matcher.group(3));
			} else {
				throw new IOException("knapsack data file not properly formatted: invalid specification line");
			}
		
			capacity = new int[nsacks];
			profit = new int[nobjectives][nitems];
			weight = new int[nsacks][nitems];
			itemNames = new String[nitems];
		
			for (int i = 0; i < nsacks; i++) {
				line = lineReader.readLine(); // line containing "="
				line = lineReader.readLine(); // line containing "knapsack i:"
				line = lineReader.readLine(); // the knapsack capacity
				matcher = capacityLine.matcher(line);
					
				if (matcher.matches()) {
					capacity[i] = Integer.parseInt(matcher.group(1));
				} else {
					throw new IOException("knapsack data file not properly formatted: invalid capacity line");
				}
		
				for (int j = 0; j < nitems; j++) {
					line = lineReader.readLine(); // line containing "item j:"
					line = lineReader.readLine(); // the item weight
					matcher = weightLine.matcher(line);
						
					if (matcher.matches()) {
						weight[i][j] = Integer.parseInt(matcher.group(1));
					} else {
						throw new IOException("knapsack data file not properly formatted: invalid weight line");
					}
								
				}
			}
				
				for (int i= 0; i < nobjectives; i++ ){
					line = lineReader.readLine(); // line containing "="
					line = lineReader.readLine(); // line containing "objective i:"
				
					for (int j = 0; j < nitems; j++) {
						line = lineReader.readLine(); // line containing "item j:"
						line = lineReader.readLine(); // the item profit
						matcher = profitLine.matcher(line);
						
						if (matcher.matches()) {
							profit[i][j] = Integer.parseInt(matcher.group(1));
						} else {
							throw new IOException("knapsack data file not properly formatted: invalid profit line");
						}
					}		
				
				}
				line = lineReader.readLine(); // line containing "Item names:"
				for (int i= 0; i < nitems; i++ ){				
					line = lineReader.readLine(); // line containing i'th items name
					itemNames[i] = line;
				}
			} finally {
				if (lineReader != null) {
					lineReader.close();
				}
			}
	}

	/**
	 * Starts the example running the knapsack problem.
	 * 
	 * @param args the command line arguments
	 * @throws IOException if an I/O error occurred
	 */
	public void printKnapsacks(Solution solution){
		boolean[] d = EncodingUtils.getBinary(solution.getVariable(0));
		double[] f = new double[nobjectives];
		double[] g = new double[nsacks];

		// calculate the  weights for the knapsacks
		for (int i = 0; i < nitems; i++) {
			if (d[i]) {
				for (int j = 0; j < nsacks; j++) {
					g[j] += weight[j][i];
				}
			}
		}
		
		for (int i = 0 ; i < nsacks; i++){
			System.out.println("Knapsack "+ (i+1)+ "= "+g[i]+"/"+capacity[i]);
		}
	}
	
	public void printMenu (Solution solution){
		boolean[] d = EncodingUtils.getBinary(solution.getVariable(0));
		System.out.print("Menu: ");
		for (int k = 0 ; k<nitems; k++){
			if (d[k])
				System.out.print("("+(k+1)+ ")+");
		}
		System.out.println();
		System.out.print("Menu: ");
		// calculate the  weights for the knapsacks
		for (int i = 0; i < nitems; i++) {
			if (d[i]) {
				System.out.print("("+itemNames[i]+")+");
			}
		}
		System.out.println();
	}
	
/*	public static void main(String[] args) throws IOException {

	String fileName = "knapsack.283.11.3"; //standard problem'daki fileName2 de deðiþtir
	String path = "src/org/moeaframework/problem/knapsack/"+ fileName;
	 NondominatedPopulation result = new Executor()
			 .withProblemClass(Knapsack.class,
					 new File(path))
			 .withAlgorithm("NSGA2")
			 .withMaxEvaluations(100000)
			 .distributeOnAllCores()
			 .run();
	
	 System.out.println("Eventually we have " + result.size() +" solutions!");
	 KnapsackExample knp = new KnapsackExample(path);
	 for (int i = 0; i < result.size(); i++) {
		 Solution solution = result.get(i);
		 System.out.println("-----Solution " + (i+1) + ":"+ solution.getVariable(0));
		 double[] objectives = solution.getObjectives();
		
		 // negate objectives to return them to their maximized
		 // form
		 objectives = org.moeaframework.util.Vector.negate(objectives);
		 
		
		 for (int j = 0 ; j < solution.getNumberOfObjectives(); j++){
			 System.out.println(" Objective "+(j+1)+" Profit: " + objectives[j]);;		 
		 }
		 if (solution.violatesConstraints())
			 System.out.println("Infeasible solution :(");
		 else 
			 System.out.println("Feasible solution :)");
		 knp.printKnapsacks(solution);
		 knp.printMenu(solution);
	
	 }
	 System.out.println("---------------------------");		 
	 

*/
	public static void main(String[] args) throws IOException {
			
		String[] algorithms = { "NSGAII", "NSGA3", "SPEA2" };
		int run_count = 10; //minimum 20 run (seed)
		int evaluation_count = 100000; //minimum 10.000 evaluations
		String fileName = "knapsack.283.11.3"; //standard problem'daki fileName2 de deðiþtir
		
		String path = "src/org/moeaframework/problem/Knapsack/"+ fileName;
		
		String firstRef = "NSGA2";
		String secondRef = "SPEA2";
		String thirdRef = "NSGA3";
	//	String fourthRef = "SMS-EMOA";
		
		
		NondominatedPopulation firstResults = new NondominatedPopulation(); 
		NondominatedPopulation secondResults = new NondominatedPopulation(); 
		NondominatedPopulation thirdResults = new NondominatedPopulation(); 
		//NondominatedPopulation fourthResults = new NondominatedPopulation(); 
		ReferenceSetMerger merger = new ReferenceSetMerger();
		
		firstResults = new Executor()
		                .withProblemClass(Knapsack.class, new File (path))
		                .withAlgorithm(firstRef)
		                .withMaxEvaluations(evaluation_count)
		                .run();
			
		merger = new ReferenceSetMerger();
		merger.add(firstRef, firstResults);
		
		System.out.println(firstRef + " is merged");
		secondResults = new Executor()
		                .withProblemClass(Knapsack.class, new File (path))
		                .withAlgorithm(secondRef)
		                .withMaxEvaluations(evaluation_count)
		                .run();		
		merger = new ReferenceSetMerger();
		merger.add(secondRef, secondResults);
		System.out.println(secondRef + " is merged");
		
		thirdResults = new Executor()
                .withProblemClass(Knapsack.class, new File (path))
                .withAlgorithm(thirdRef)
                .withMaxEvaluations(evaluation_count)
                .run();		
		merger = new ReferenceSetMerger();
		merger.add(thirdRef, thirdResults);
		System.out.println(thirdRef + " is merged");
		

		
		
	    new Plot()
         	.add(firstRef, firstResults)
            .show();
	   
	    PopulationIO.writeObjectives(new File(fileName + ".pf"), merger.getCombinedPopulation());
   	   
		//setup the experiment
		Executor executor = new Executor()
			.withProblemClass(Knapsack.class, new File (path))
		    .withMaxEvaluations(evaluation_count);
	
		Analyzer analyzer = new Analyzer().withSameProblemAs(executor)
			.withSameProblemAs(executor)
			.includeHypervolume()
		    .showStatisticalSignificance();
		
		//run each algorithm for run_count times seeds
		for (String algorithm : algorithms) {
		    analyzer.addAll(algorithm, executor.withAlgorithm(algorithm).runSeeds(run_count));  
			System.out.println(algorithm+ " is run!");
		}
		
		analyzer.printAnalysis();

	} 
}

	