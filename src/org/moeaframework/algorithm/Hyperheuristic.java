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
package org.moeaframework.algorithm;



import java.io.File;

import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.algorithm.GDE3;
import org.moeaframework.algorithm.NSGAII;
import org.moeaframework.algorithm.single.AggregateObjectiveComparator;
import org.moeaframework.algorithm.single.LinearDominanceComparator;
import org.moeaframework.core.Algorithm;
import org.moeaframework.core.Initialization;

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
import org.moeaframework.core.operator.real.DifferentialEvolution;
import org.moeaframework.core.operator.real.DifferentialEvolutionSelection;
import org.moeaframework.core.operator.real.PM;
import org.moeaframework.core.operator.real.SBX;
import org.moeaframework.core.spi.OperatorFactory;
import org.moeaframework.problem.chemicalShiftAssignment.RunHyperheuristic;
import org.moeaframework.problem.chemicalShiftAssignment.chemicalShiftAssignment;
import org.moeaframework.util.TypedProperties;

@SuppressWarnings("deprecation")
/*
 * Public class to define the hyperheuristic algorithm. 
 */
public class Hyperheuristic extends AbstractEvolutionaryAlgorithm {

	private NSGAII nsga2;
	private NSGAII nsga3;
	private Algorithm ga;
	private HillClimber hc;
	
	private int iteration;

	/*
	 * The algorithms that will be used in the hyperheuristic algorihm are defined. 
	 */
	public Hyperheuristic(chemicalShiftAssignment problem, NondominatedSortingPopulation population, Initialization initialization,double crossoverProbability, double mutationProbability) {
		
		super(problem, population, null, initialization);
		
		
		//Hill Climber
		Variation variationHillClimber = new PM(1.0 / problem.getNumberOfVariables(),30.0);
		hc = new HillClimber(problem, population, initialization, variationHillClimber);
		
		//NSGA2
		Variation variation = new GAVariation(new SBX(crossoverProbability, 25.0),new PM(mutationProbability / problem.getNumberOfVariables(), 30.0));
		TournamentSelection selectionNSGA2 = new TournamentSelection(2,new ChainedComparator(new ParetoDominanceComparator(),new CrowdingComparator()));
		nsga2 = new NSGAII(problem, population, null, selectionNSGA2, variation, initialization); //null = no archive

		//GA
		AggregateObjectiveComparator comparator = new LinearDominanceComparator(new double[] { 1.0 }); // default value taken from StandardAlgorithms.f90 newGeneticAlgorithm
		TournamentSelection selectionGA = new TournamentSelection(2, comparator);
		ga = new org.moeaframework.algorithm.single.GeneticAlgorithm(problem, comparator, initialization, selectionGA, variation);
		
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
		population = new ReferencePointNondominatedSortingPopulation(problem.getNumberOfObjectives(), divisionsOuter, divisionsInner);
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
		nsga3 =  new NSGAII(problem, population, null, selectionNSGA3, variationNSGA3, initialization);
		//NSGA3 end

	}

	/*
	 * The iterate class for the hyperheuristic algorithm.
	 * 
	 */
	@Override
	protected void iterate() {
		
	
		if (iteration % 2 == 0) {
			nsga2.iterate();
		} else {
	
			hc.iterate();
			
		}
		numberOfEvaluations = nsga2.getNumberOfEvaluations() +	hc.getNumberOfEvaluations();
	
		iteration++;
		
	} 
	
	/*
	 * This function prints the given population
	 */
	void printPopulation (NondominatedSortingPopulation population) {
		int i, j;
		
		System.out.println("----------------POPULATION------------------------");
		for (i = 0; i < population.size(); i++) {
			System.out.print(i+1+". individual: ");
			for (j = 0; j < population.get(i).getNumberOfVariables(); j++) {
				System.out.print(population.get(i).getVariable(j)+ " ");
			}
			System.out.println();
		}
		
		System.out.println("--------------------------------------------------");
		
	}

}
