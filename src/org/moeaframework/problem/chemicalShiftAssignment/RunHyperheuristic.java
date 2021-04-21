/* Copyright 2016-2021 EMY
 *
 * This file is a part of the Chemical Shift Assignment Solution for EMY. 
 *
 * These classes are added into the MOEA Framework, where the algorithms are implemented.
 * 
 * The MOEA Framework and this solution software is free software:you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or (at your
 * option) any later version.
 *
 * 
 */

package org.moeaframework.problem.chemicalShiftAssignment;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.nio.file.StandardOpenOption;

import org.moeaframework.algorithm.HillClimber;
import org.moeaframework.algorithm.Hyperheuristic;
import org.moeaframework.algorithm.NSGAII;
import org.moeaframework.algorithm.ReferencePointNondominatedSortingPopulation;
import org.moeaframework.algorithm.single.AggregateObjectiveComparator;
import org.moeaframework.algorithm.single.LinearDominanceComparator;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedPopulation;
import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Selection;
import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.core.comparator.AggregateConstraintComparator;
import org.moeaframework.core.comparator.ChainedComparator;
import org.moeaframework.core.comparator.CrowdingComparator;
import org.moeaframework.core.comparator.DominanceComparator;
import org.moeaframework.core.comparator.ParetoDominanceComparator;
import org.moeaframework.core.operator.GAVariation;
import org.moeaframework.core.operator.RandomInitialization;
import org.moeaframework.core.operator.TournamentSelection;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.util.TypedProperties;

/*
 * Public class to run the main function that defines the configuration parameters below for the algorithms:
 * 	- runCount: How many times one algorithm to be run. The default value is set to "1".
 * 	- tuningCount: How many times tuning (=repetition parameter) should be done in one run of the algorithm. The default value is set to "10".
 * 	- evaluationCount: How many times one algorithm should be iterated. This is basically the iteration count for evolutionary algorithm. The default value is set to "10.000".
 * 	- algorithmName: Which algorithm to be run. Options: NSGA2, NSGA3, GA or Hyperheuristic. The default value is set to "GA".
 * 	- populationSize: The size of the population for evolutionary optimization algorithms. The default value is set to "50".
 *	- crossoverProbability: The probability of applying crossover in each iteration. The default value is set to "1.0".
 *	- mutationProbability: The probability of applying mutation in each iteration. The default value is set to "1.0".
 * 	- initMode: The method of initializing the population for the evolutionary optimization algorith. Options: manual, random. The default value is "random".
 * 
 */
public class RunHyperheuristic {
	private String fileName  = ""; 
	private static String p = "src/org/moeaframework/problem/chemicalShiftAssignment/configuration.txt"; //all input file names are in this file
	private String path = p + fileName;
	private static String rootResultPath = "ChemicalShiftAssignmentResults/"; //IMPORTANT! Files to be deleted, so do not change this path.
		
			
	public static void main(String[] args) throws IOException {
		int i, j, k, r;
		int randomCounter = 0;
		double best = -Double.MAX_VALUE;
		
		/*
		 * Configuration Values
		 */
		int runCount = 1;
		int tuningCount = 10; 
		int evaluationCount = 10050; // First 50 evaluations are used to initialize (if #pop is 50, then 50 evaluations for that)
		String algorithmName = "GA"; //GA NSGA2 NSGA3 Hyperheuristic
		int populationSize = 50; //50;  //number of individuals in one population
		double crossoverProbability = 1.0;
		double mutationProbability = 1.0;
		String initMode = "manual"; // manual random Default: random
		boolean writeObjectivesToFile = false;
		
		/*
		 * END 
		 */

		long startTime, startTimeOverall, durationOverall;
		RunHyperheuristic hh = new RunHyperheuristic();
		File directory = new File (rootResultPath);
		emptyResultDirectory(directory);
		
		
		System.out.println("Source code of Emel for GitHub...");
		for (r = 0; r < runCount; r++){
			
			System.out.println("********************RUN ("+ (r +1) +")*************************************");
			String resultPath = rootResultPath + "/"+r+".Run/";
			
			startTimeOverall = System.nanoTime();
			for (j = 0; j < tuningCount; j++){
				startTime = System.nanoTime();
				System.out.println ("==================TUNING ("+ (j +1) +")===("+(r +1)+"/"+runCount+" run)=====================");
				System.out.println("Manual Run with "+ algorithmName +" and "+evaluationCount +" evaluations.");
					
				// define the problem
				chemicalShiftAssignment problem = new chemicalShiftAssignment(new File(hh.getPath()));
	
				Initialization initialization;
				if (initMode.equals("manual")) { //manual initialization
					initialization = new ChemicalShiftAssignmentInitialization(problem,populationSize);
				} else { //random (default)
					initialization = new RandomInitialization(problem, populationSize);
				}
				
				
				
				// define the selection operator
				TournamentSelection selection = new TournamentSelection(2,new ChainedComparator(new ParetoDominanceComparator(),new CrowdingComparator()));
					
				// define the crossover / mutation operator
				Variation variation = new GAVariation(new SBX(crossoverProbability, 25.0),new PM(mutationProbability / problem.getNumberOfVariables(), 30.0));
	
				// construct the algorithm
				//NSGA2
				Algorithm nsga2 = new NSGAII(problem,new NondominatedSortingPopulation(), null, selection, variation, initialization); //null is for no archive
				
				//GA
				AggregateObjectiveComparator comparator = new LinearDominanceComparator(new double[] { 1.0 }); // default value taken from StandardAlgorithms.f90 newGeneticAlgorithm
				selection = new TournamentSelection(2, comparator);
				Algorithm ga = new org.moeaframework.algorithm.single.GeneticAlgorithm(problem, comparator, initialization, selection, variation);
	
				//NSGA3
				int divisionsOuter = 4;
				int divisionsInner = 0;
				TypedProperties properties = new TypedProperties();				
				if (problem.getNumberOfObjectives() == 1) {
					divisionsOuter = 100;
				} else if (problem.getNumberOfObjectives() == 2) {
					divisionsOuter = 99;
				} else if (problem.getNumberOfObjectives() == 3) {
					divisionsOuter = 12;
				} else if (problem.getNumberOfObjectives() == 4) {
					divisionsOuter = 8;
				} else if (problem.getNumberOfObjectives() == 5) {
					divisionsOuter = 6;
				} else if (problem.getNumberOfObjectives() == 6) {
					divisionsOuter = 4;
					divisionsInner = 1;
				} else if (problem.getNumberOfObjectives() == 7) {
					divisionsOuter = 3;
					divisionsInner = 2;
				} else if (problem.getNumberOfObjectives() == 8) {
					divisionsOuter = 3;
					divisionsInner = 2;
				} else if (problem.getNumberOfObjectives() == 9) {
					divisionsOuter = 3;
					divisionsInner = 2;
				} else if (problem.getNumberOfObjectives() == 10) {
					divisionsOuter = 3;
					divisionsInner = 2;
				} else {
					divisionsOuter = 2;
					divisionsInner = 1;
				}
				ReferencePointNondominatedSortingPopulation population = new ReferencePointNondominatedSortingPopulation(problem.getNumberOfObjectives(), divisionsOuter, divisionsInner);
				Selection selectionNSGA3 = null;
				if (problem.getNumberOfConstraints() == 0) {
					selectionNSGA3 = new Selection() {		
						@Override
						public Solution[] select(int arity, org.moeaframework.core.Population population) {
							Solution[] result = new Solution[arity];							
							for (int i = 0; i < arity; i++) {
								result[i] = population.get(PRNG.nextInt(population.size()));
							}							
							return result;
						}				
					};
				} else {
					selectionNSGA3 = new TournamentSelection(2, new ChainedComparator(
							new AggregateConstraintComparator(),
							new DominanceComparator() {
								@Override
								public int compare(Solution solution1, Solution solution2) {
									return PRNG.nextBoolean() ? -1 : 1;
								}							
							}));
				}
				properties.setBoolean("sbx.swap", false);
				properties.setDouble("sbx.distributionIndex", 30.0);
				properties.setDouble("pm.distributionIndex", 20.0);	
				properties.setDouble("hux.rate", crossoverProbability);
				properties.setDouble("bf.rate", mutationProbability);
				Variation variationNSGA3 = OperatorFactory.getInstance().getVariation(null, properties, problem);
				Algorithm nsga3 =  new NSGAII(problem, population, null, selectionNSGA3, variationNSGA3, initialization);
				//NSGA3 end
	
				// run the algorithm for 10,000 evaluations	
				if (algorithmName == "NSGA2"){	
					System.out.println(">>RunHyperheuristic: Run NSGA2 algorithm with "+evaluationCount+ " evaluation and "+tuningCount+" runs!");
					while (nsga2.getNumberOfEvaluations() < evaluationCount) {
						nsga2.step();
						writeObjectivesToFile(nsga2, writeObjectivesToFile, nsga2.getNumberOfEvaluations());						
					}
					best = evaluateResult( nsga2.getResult(), r, j, best, resultPath);
					
				} else if (algorithmName == "GA"){	
					System.out.println(">>RunHyperheuristic: Run GA algorithm with "+evaluationCount+ " evaluation and "+tuningCount+" runs!");	
					while (ga.getNumberOfEvaluations() < evaluationCount) {
						ga.step();
						writeObjectivesToFile(ga, writeObjectivesToFile, ga.getNumberOfEvaluations());
						
					}
					best = evaluateResult(ga.getResult(), r, j, best, resultPath);
					
				} else if (algorithmName == "NSGA3"){
					System.out.println(">>RunHyperheuristic: Run NSGA3 algorithm with "+evaluationCount+ " evaluation and "+tuningCount+" runs!");	
					while (nsga3.getNumberOfEvaluations() < evaluationCount) {
						nsga3.step();
						writeObjectivesToFile(nsga3, writeObjectivesToFile, nsga3.getNumberOfEvaluations());
					}	
					best = evaluateResult(nsga3.getResult(), r, j, best, resultPath);
				} else if (algorithmName == "Hyperheuristic") {
					System.out.println(">>RunHyperheuristic: Run Hyperheuristic algorithm with "+evaluationCount+ " evaluation and "+tuningCount+" runs!");	
					Algorithm hyperheuristic = new Hyperheuristic (problem, new NondominatedSortingPopulation(), initialization, crossoverProbability, mutationProbability);
					System.out.println("Defined evaluationCount: " + evaluationCount);
					System.out.println("Current evaluationCount: "+ hyperheuristic.getNumberOfEvaluations());
					while (hyperheuristic.getNumberOfEvaluations() < evaluationCount) {
						hyperheuristic.step();
						System.out.println("after one step: current evaluation count: "+ hyperheuristic.getNumberOfEvaluations());
					}
					best = evaluateResult(hyperheuristic.getResult(), r, j, best, resultPath);
				}
				
				final long duration = System.nanoTime() - startTime;
				System.out.println("Duration of the current tuning: "+ duration/1000000000 +" seconds");
			}
			durationOverall = (System.nanoTime() - startTimeOverall)/1000000000;
			System.out.println("Duration of all runs: "+ durationOverall +" seconds");
	
			System.out.println("****Best score: "+ best);
			writeDurationsToFile (rootResultPath, r, durationOverall);
		}
	}
	
	/*
	 * This function writes the objective values for each of them to the given file.
	 * This function is called only if it is explicitly turned on in the configuration part in the main function.
	 */
	private static void writeObjectivesToFile(Algorithm alg, boolean writeToFile, int evaluationCount){
		int i, j, popSize, numberOfObjectives;
		if (!writeToFile)
			return;
		popSize = alg.getResult().size();
		numberOfObjectives = alg.getProblem().getNumberOfObjectives();
		System.out.println("===================="+ evaluationCount+". evaluation===================================");
		for (i = 0; i < popSize; i++) {
			for (j = 0 ; j < numberOfObjectives; j++) {
				System.out.print(alg.getResult().get(i).getObjective(j)+ "\t");
			}
			System.out.println();
			
			
		}
	}
	
	/*
	 * This function writes the duration of each run to the predefined text file.
	 */
	private static void writeDurationsToFile (String rootResultPath, int r, final long durationOverall) throws IOException{

		String fileName = rootResultPath +"Durations.txt";
		
		try{	
			File file = new File(fileName);
			FileWriter fr = new FileWriter(file, true);
			BufferedWriter br = new BufferedWriter(fr);
			PrintWriter writer = new PrintWriter(br);
			writer.println(String.format("%d.RUN: %-10s seconds.", (r+1), durationOverall));
			
			writer.close();
			
		}  finally {
			
			
		}
	
		
	}
	
	/*
	 * This function evaluates the final results of the objectives. 
	 */
	private static double evaluateResult(NondominatedPopulation result, int r, int j, double best, String resultPath) throws IOException {
		
		int i, k;
		double G;
		
		System.out.println("Eventually we have " + result.size() +" solutions!");
	
		int[] asg = new int[result.get(0).getNumberOfVariables()];
		
		for (i = 0; i < result.size(); i++) {  //check all results one by one
			 Solution solution = result.get(i);
			 System.out.println("-----Solution " + (i+1) + "-----------------");
			 for (k = 0; k < solution.getNumberOfVariables(); k++){
				 System.out.print(solution.getVariable(k)+" ");
				 asg[k] = Integer.parseInt(solution.getVariable(k).toString());
			 }
			 System.out.println();
			 G = solution.getG();
			 System.out.println("G= "+ G);

			 double[] objectives = solution.getObjectives();
			
			 //negate objectives to return them into their maximized form
			 //objectives are: -a, -b, d and -c. So, only the first two and last is negated here:
			 objectives[0] = -objectives[0]; 
			 objectives[1] = -objectives[1];
			 objectives[3] = -objectives[3];
			 System.out.println("a: "+ objectives[0]);
			 System.out.println("b: "+ objectives[1]);
			 System.out.println("d: "+ objectives[2]);
			 System.out.println("c: "+ objectives[3]);

			 solution.writeResultsToFile(rootResultPath, r, j, i);
			 
		
			//shows the best result only for the first objective, others value is ignored
			 if (G > best){
					 best = G;
					 System.out.println("best assigned to: "+ best);
			}
	
			 if (solution.violatesConstraints())
				 System.out.println("Infeasible :(");
			 else 
				 System.out.println("Feasible :)");
		 } 
		
		return best;
		
	}

	/*
	 * This function deletes the existing results from previous runs once the algorithm starts. 
	 */
	private static void emptyResultDirectory(File directory)  {
		File[] files;
		File currentFile;
		int i; 
		
		if (directory.isDirectory()) {
		    files = directory.listFiles();
		    for (i = 0 ; i <  files.length; i++){
		    	currentFile = files[i];
		    	emptyResultDirectory(currentFile);
		    }

		  }
		  try{
			 directory.delete();		
		  } finally{
			  
		  }

	}

	public String getPath() {
		
		return path;
	}

	public void setPath(String path) {
	
		this.path = path;
	}
}