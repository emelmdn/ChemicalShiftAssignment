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
 * Public class Tree to save and operate on the search space with building a tree structure of the nodes.
 */
public class Tree {
	private int measuredPeakNumber; //number of measured peak of the spectrum
	private ArrayList<Integer> orientation;
	private ArrayList<Double> measuredPeakInSpectrum; //number of measured peak of the spectrum
	private ArrayList<Tree> followingDimension; //number of measured peak of the spectrum
	private ArrayList<Tree> furtherDimension; //number of measured peak of the spectrum

	
	public Tree() {
	
		this.measuredPeakNumber = 0;
		this.orientation = new ArrayList<Integer>();
		this.measuredPeakInSpectrum = new ArrayList<Double>();
		this.followingDimension = new ArrayList<Tree>();
		this.furtherDimension = new ArrayList<Tree>();
	}

	/*
	 * copy constructor
	 * 
	 */
	public Tree(Tree r) {
		int i,j;
		
		this.measuredPeakNumber = r.measuredPeakNumber;
		
		this.orientation = new ArrayList<Integer>(r.orientation.size());
		for (i = 0; i < r.orientation.size(); i++)
			this.orientation.add(r.orientation.get(i));
		
		this.measuredPeakInSpectrum = new ArrayList<Double>(r.measuredPeakInSpectrum.size());
		for (i = 0; i < r.measuredPeakInSpectrum.size(); i++)
			this.measuredPeakInSpectrum.add(r.measuredPeakInSpectrum.get(i));
		
		this.followingDimension = new ArrayList<Tree>(r.followingDimension.size());	
		for (i = 0; i < r.followingDimension.size(); i++){
			this.followingDimension.add(r.followingDimension.get(i));
		}
		
		this.furtherDimension = new ArrayList<Tree>(r.furtherDimension.size());
		for (i = 0; i < r.furtherDimension.size(); i++){
			this.furtherDimension.add(r.furtherDimension.get(i));
		}
	
	}
	

	
	public ArrayList<Integer> getOrientation(){
		return this.orientation;
	}
	
	public int getOrientationOf(int i){
		return this.orientation.get(i);
	}
	
	public double getMeasuredPeakInSpectrum(int i){
		return this.measuredPeakInSpectrum.get(i);
	}


	public ArrayList<Double>  getMeasuredPeaks(){
		return this.measuredPeakInSpectrum;
	}
	
	
	
	public void setDimensions(double[] array){
		int i;
		
		int size = array.length;
		
		this.orientation = new ArrayList<Integer>(1);
		this.measuredPeakInSpectrum.clear();
		this.measuredPeakInSpectrum = new ArrayList<Double>(size);
		for (i = 0; i < size; i++){
			this.measuredPeakInSpectrum.add(array[i]);
		}
		this.followingDimension = new ArrayList<Tree>(size);
		for (i = 0; i < size; i++)
			this.followingDimension.add(new Tree());
		
		this.furtherDimension = new ArrayList<Tree>(size);
		for (i = 0; i < size; i++)
			this.furtherDimension.add(new Tree());
	}
	
	public void setFollowingDimensionToSingel(double value){
		this.measuredPeakInSpectrum.clear();
		this.measuredPeakInSpectrum.add(value);
		this.followingDimension.add(new Tree());
		this.furtherDimension.add(new Tree());
		this.measuredPeakNumber = 1;
		
	}
		
	public Tree getFollowingDimension(int i){
		return this.followingDimension.get(i);
	}
	
	public Tree getFurtherDimension(int i){
		return this.furtherDimension.get(i);
	}
	public void setMeasuredPeakNumber(int i){
		this.measuredPeakNumber = i;
	}
	public int getMeasuredPeakNumber(){
		return this.measuredPeakNumber;
	}

	public void setOrientationTo(int[] input) {
		int i;
		
		this.orientation = new ArrayList<Integer>(input.length);
		for (i = 0; i < input.length; i++){
			this.orientation.add(input[i]);
		}
		
	}
	public void addToOrientation(int input) {
		orientation.add(input);
		
	}


}
