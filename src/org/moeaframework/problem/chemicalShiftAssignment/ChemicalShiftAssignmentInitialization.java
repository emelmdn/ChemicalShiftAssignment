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

import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.moeaframework.util.io.CommentedLineReader;

import org.moeaframework.core.Initialization;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryIntegerVariable;
import org.moeaframework.util.io.CommentedLineReader;
import org.moeaframework.util.sequence.LatinHypercube;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.problem.chemicalShiftAssignment.*;

/*
 * Public class to implement the initialization method for the population of evolutionary optimization algorithm.
 */
public class ChemicalShiftAssignmentInitialization implements Initialization{
	
	String path = "src/org/moeaframework/problem/chemicalShiftAssignment/";
	private chemicalShiftAssignment problem;
	private int size;  //number of individuals in one population
	private int variableCount; //number of variables in one individual(=#expected peaks)
	
	public ChemicalShiftAssignmentInitialization(chemicalShiftAssignment problem, int populationSize) {
		super();
		this.problem = problem;
		this.problem.setPopulationSize(populationSize);
		this.size = populationSize;
		this.variableCount = problem.getNumberOfVariables();

	}
	
	/*
	 * This function initializes all individuals manually
	 *
	 */
	public Solution[] initialize(){
		
		
				int k;
				Solution[] result = new Solution[size];  //one population, including "size" number of individuals
				final long startTime = System.nanoTime(); 
				System.out.println("ALL INDIVIDUALS WILL BE CREATED MANUALLY!");
				for (int i = 0; i < size; i++) {  //create "size" number of individuals in one population
						
					Solution solution = problem.newSolution();
					solution = problem.initializeOneIndividual(i);
					
					int[] asg = new int[solution.getNumberOfVariables()];
					
					 for (k = 0; k < solution.getNumberOfVariables(); k++){
						 asg[k] = Integer.parseInt(solution.getVariable(k).toString());
					
					 }
					 
					result[i] = solution;
	
				}
				final long currentDuration = System.nanoTime() - startTime;
				System.out.println("Duration until individuals are created: "+ currentDuration/1000000000 +" seconds");
				return result;
		}
				 
		

	/*
	 * This function creates half of the population manually, the rest randomly.
	 */
	public Solution[] initialize_halfhalf() {
		int k, i, j;
		int half;
		
		half = size / 2;
		
	
		Solution[] result = new Solution[size];  //one population, including "size" number of individuals
		System.out.println("HALF OF THE INDIVIDUALS WILL BE CREATED RANDOMLY, OTHER HALF MANUALLY!");
		for (i = 0; i < half; i++) {  //create "size" number of individuals in one population
				
			Solution solution = problem.newSolution();
			solution = problem.initializeOneIndividual(i);
			
			
			int[] asg = new int[solution.getNumberOfVariables()];
			
			 for (k = 0; k < solution.getNumberOfVariables(); k++){
				 System.out.print(solution.getVariable(k)+" ");
				 asg[k] = Integer.parseInt(solution.getVariable(k).toString());
			
			 }
				 
	
			result[i] = solution;
			
		}
		for (i = half; i < size; i++) {  //create "size" number of individuals in one population
				
			Solution solution = problem.newSolution();
			
			for (j = 0; j < solution.getNumberOfVariables(); j++) {
				solution.getVariable(j).randomize();
			}
			
			
			result[i] = solution;

			int[] asg = new int[solution.getNumberOfVariables()];
			
			 for (k = 0; k < solution.getNumberOfVariables(); k++){
				 System.out.print(solution.getVariable(k)+" ");
				 asg[k] = Integer.parseInt(solution.getVariable(k).toString());
			
			 }
		
			result[i] = solution;
			
		}
		return result;
		
		 
	}
	
	
}

 