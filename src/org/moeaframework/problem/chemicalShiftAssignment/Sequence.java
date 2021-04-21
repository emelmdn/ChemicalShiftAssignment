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

import java.util.ArrayList;

/*
 * Public class to save the Sequence values of the NMR experiments
 */
public class Sequence {
	int residueNameLength = 5;
	int number;
    int[] residueIds;
    String[] externalNumber = new String[residueNameLength];
    double[] score;
    double[] scoreResidue;
    double[] averageScore;
    double[] deviation;
    String[] namesResidue = new String[residueNameLength];
    boolean averageIsSet;
    double[] connectSequence;
   
	public Sequence(ArrayList<String> aminoacidSequence) {
		int i;
		
		this.number = aminoacidSequence.size();
		score = new double[number];
		scoreResidue = new double[number];  
		averageScore = new double[number];
        deviation = new double[number];
        
        residueIds = new int[number];
        
        namesResidue = new String[number];
        connectSequence = new double[number];
        externalNumber = new String[number];
		
        score = new double[number];
		score = new double[number];
		score = new double[number];
		score = new double[number];
		score = new double[number];
        
		for (i = 0; i < connectSequence.length; i++)
			this.connectSequence[i] = 0;
  		for (i = 0 ; i < averageScore.length; i++) {
			averageScore[i] = 1.0;
		}
		for (i = 0; i < deviation.length; i++)
			deviation[i] = 0.2;
		averageIsSet = false;
    }

	
	public double getScore(int index) {
		return this.score[index];
	}
	
	public double getConnectSequence(int index) {
		return this.connectSequence[index];
	}
	public boolean getAverageIsSet() {
		return this.averageIsSet;
	}
	public void setAverageIsSet(boolean newValue) {
		this.averageIsSet = newValue;
	}
	public double getAverageScore(int index) {
		return this.averageScore[index];
	}
	public void setAverageScore(int index, double value) {
		this.averageScore[index] = value;
	}
	
	public double getDeviation(int index) {
		return this.deviation[index];
	}
	
	public int getDeviationLength() {
		return this.deviation.length;
	}
	
	public void setDeviation(int index, double value) {
		this.deviation[index] = value;
	}
	
	public int getAverageScoreLength() {
		return this.averageScore.length;
	}
	public int getNumber() {
		return this.number;
	}
	public double getScoreResidue(int index) {
		return this.scoreResidue[index];
	}
}

