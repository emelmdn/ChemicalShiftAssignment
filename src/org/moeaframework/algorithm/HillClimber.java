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


import org.moeaframework.algorithm.AbstractEvolutionaryAlgorithm;
import org.moeaframework.core.Initialization;
import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.PRNG;
import org.moeaframework.core.Problem;

import org.moeaframework.core.Solution;
import org.moeaframework.core.Variation;
import org.moeaframework.problem.chemicalShiftAssignment.chemicalShiftAssignment;

/*
 * Public class that implements the standard hill climber algorithm. 
 */
public class HillClimber extends AbstractEvolutionaryAlgorithm {
	
	private final Variation variation;
	private chemicalShiftAssignment problem;
	
	public HillClimber(chemicalShiftAssignment problem, NondominatedSortingPopulation population, Initialization initialization,Variation variation) {
		
		super(problem, population, null, initialization);
		this.variation = variation;
		this.problem = problem;
	}
	
	@Override
	protected void iterate() {
		
		// get the current population
		NondominatedSortingPopulation population = (NondominatedSortingPopulation) getPopulation();
		
		// randomly select a solution from the population
		int index = PRNG.nextInt(population.size());
		Solution parent = population.get(index);
	
		
		// hill climbing on the selected solution
		Solution offspring = this.problem.runHillClimber(parent);
	
		// evaluate the objectives/constraints
		evaluate(offspring);

		// add the offspring to the population
		population.add(offspring);
	

		// use non-dominated sorting to remove the worst solution
		population.remove(index); //ben worst yerine hill climbing yaptigimi cikardim
		
		
	}
	

}