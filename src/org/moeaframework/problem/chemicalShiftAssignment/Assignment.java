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
import java.util.Random;

/*
 * Public class for the chemical shift assignment of the expected peaks onto the measured peaks. 
 */
public class Assignment{
	
	private int size; //size of the assignment
	private int[] assignment;  //getStatus method below returns// -1: unassigned, > -1: assigned (shows assigned Measured Peak Id)
	private int[] rejectNumber;  //0: not rejected, >0: rejected  //if assignment value for that expPeak is -1, then check rejectNr here
	private ArrayList<Integer> assigned;
	private ArrayList<Integer> unassigned;
	private ArrayList<Integer> rejected;
	
 	long seed = System.currentTimeMillis();
	Random r = new Random (seed);
	
	
	public Assignment(){
		unassigned = new ArrayList<Integer>();
		assigned = new ArrayList<Integer>();
		rejected = new ArrayList<Integer>();
	}
	
	public Assignment(int size){
		unassigned = new ArrayList<Integer>();
		assigned = new ArrayList<Integer>();
		rejected = new ArrayList<Integer>();
		assignment = new int[size];
		rejectNumber = new int[size];
		this.size = size;
	}
	
	/*
	 * This function returns the status of the peak.
	 * Options are PF: PeakFixed, PA: PeakAssigned, PU:PeakUnassigned, PR: PeakRejected 
	 */
	public String getStatus (int expId) {
		
		if (this.assignment[expId] > -1) {
			return "PA";  //peak assigned
		} else { //assignment value is -1
			if (this.rejectNumber[expId] > 0) 
				return "PR";  //peak rejected
			else 
				return "PU";  //peak unassigned
		}
	
	}
	
	/*
	 * This function checks the assignment and returns true if there is at least one unassigned expected peak.
	 * It returns false otherwise. 
	 */
	public boolean unassignedExists(){
	
		if (this.unassigned.size() > 0)
			return true;
		else 
			return false;
	}
	
	public int[] getAssignmentWhole() {
		return this.assignment;
	}

	public int getRejectNumber (int index){
		return this.rejectNumber[index];
	}
	
	/*
	 * This function maps the given expected peak to the given measured peak. 
	 */
	public void setAssignment(int expId, int measId){
		this.assignment[expId] = measId;
	
		if (measId != -1) { //is the expId is assigned to a real measId; then remove from unassigned list. If it is set to -1, do not remove
			this.removeFromUnassigned(expId); //since expId is assigned, remove it from unassigned list
			this.addToAssigned(expId);
		}
	}
	
	/*
	 * This function adds the given expected peak to the assigned list.
	 */
	public void addToAssigned(int expId) {
		this.assigned.add(expId);
	}
	
	/*
	 * This function gets one expected peak from the unassigned list. 
	 */
	public ArrayList <Integer> getUnassigned(){
		return this.unassigned;
	}
	
	public ArrayList <Integer> getAssigned(){
		return this.assigned;
	}
	
	public void printAssignmentClass(){
		int i;
		System.out.println("---Assignment Class values:--");
		System.out.println("size: " + this.size);
		
		System.out.print("Assignment: ");
		for (i = 0; i < this.assignment.length; i++)
			System.out.print(this.assignment[i]+ " ");
		System.out.println();
		
		System.out.print("Assigned: ");
		for (i = 0; i < this.assigned.size(); i++)
			System.out.print(this.assigned.get(i)+ " ");
		System.out.println();
		
		System.out.print("Unassigned: ");
		for (i = 0; i < this.unassigned.size(); i++)
			System.out.print(this.unassigned.get(i)+ " ");
		System.out.println();
		
		System.out.println("Rejected: ");
		for (i = 0; i < this.rejected.size(); i++)
			System.out.print(this.rejected.get(i)+ " ");
		System.out.println();
		
		System.out.println("RejectNr:");
		for (i = 0; i < this.rejectNumber.length; i++)
			System.out.print(this.rejectNumber[i] + " ");
		System.out.println();
		System.out.println("-----");
	}
	
	public void printAssignment(){
		int i;
		
		System.out.print("Assignment: ");
		for (i = 0; i < this.size; i++)
			System.out.print(this.assignment[i]+ " ");
		System.out.println();
	}
	
	public void addRejectedPeak(int rejectedPeakId){
		this.rejected.add(rejectedPeakId);
		this.rejectNumber[rejectedPeakId] += 1;
	}
	
	/*
	 * This function resets all assignment values. 
	 */
	public void reset(){
		int i;
	
		this.size = this.assignment.length;
		this.unassigned = new ArrayList<Integer>(this.size);
		this.assigned.clear();
		for (i = 0; i < size; i++){
			this.assignment[i] = -1;
			this.rejectNumber[i] = 0;
			unassigned.add(i);
		}
		unassigned = listValuesRandomize(this.unassigned);
	}
	
	/*
	 * This function removes the first value from assigned list and returns it as a result
	 */
	public int dequeueFromAssignedList() {
		if (this.assigned.size() > 0)
			return this.assigned.remove(0);
		else
			return -1;
	}

	/*
	 * This function does not really dequeue the value from the list, but sends the first value and keeps it still in the list
	 */
	public int dequeueFromUnassignedList() {
		return dequeue(this.unassigned);
	}
	
	/*
	 * This function removes the first element of the list and returns it as a result
	 */
	public int dequeue(ArrayList<Integer> list){
			
			int id = -1;
			id = list.get(0);
			
			
			return id;

	}
	
	private void printList (ArrayList<Integer> asg){
			int i;
			
			for (i = 0; i < asg.size(); i++)
				System.out.print(asg.get(i) + " ");
			System.out.println();
		
	}

	/*
	 * This function randomizes the input list and returns the randomized list as output. 
	 */
	public ArrayList <Integer> listValuesRandomize (ArrayList <Integer> list){

		int beginIndex, endIndex;
		int length, positionForChange, temp, i;
		double random;
		
		ArrayList<Integer> result = new ArrayList<Integer>(list.size());
		
		int []array = new int[list.size()];

		for (i = 0; i<list.size(); i++)
			array[i] = list.get(i);
		
		beginIndex = 0;
		endIndex = list.size();
		
		length = endIndex - beginIndex; 
		
		if (length > 1) {
			for (i = beginIndex; i < endIndex; i++){
				random = this.r.nextDouble();
				positionForChange = (int) Math.floor(random * length) + beginIndex;
				temp = array[i];
				array[i] = array[positionForChange];
				array[positionForChange] = temp;		
			}
		}
		
		for (i = 0; i<array.length; i++){
			result.add(array[i]);
		}
		
		
		return result;
		
	}
	
	/*
	 * This function removes the received expected peak Id from the unassigned list, since it is assigned. 
	 */
	public void removeFromUnassigned(int expectedPeakId) {
		
		if (this.unassignedExists()){
			if (this.unassigned.get(0) == expectedPeakId)
				this.unassigned.remove(0);
			else{
				for (int i = 0; i < this.unassigned.size(); i++){
					if (this.unassigned.get(i) == expectedPeakId)
						this.unassigned.remove(i);
				}
			}
		}
		
	}

	public void resetUnassigned() {
		this.unassigned.clear();
		
	}

	public void addToUnassigned(int pExpPeak) {

		this.unassigned.add(pExpPeak);
		this.rejectNumber[pExpPeak] = 0;
		
	}

	public int dequeueRejected() {
		if (this.rejected.size() > 0)
			return this.rejected.remove(0);
		else
			return -1;
	}

	public int getAssignment(int expId) {
		return this.assignment[expId];
	}

	public void setExpStatusToUnassigned(int comExpId) {
		this.rejectNumber[comExpId] = 0;
	
	}

	public void replaceUnassignedWith(int[] vals) {
		int i;
		
		this.unassigned.clear();
		
		for (i = 0; i < vals.length; i++)
			this.unassigned.add(vals[i]);
		
	}

	public void printUnassigned() {
		int i;
		System.out.print("Unassigned: [");
		for (i = 0; i < this.unassigned.size(); i++)
			System.out.print(this.unassigned.get(i)+ " ");
		System.out.println("]");
	}
	
	public void printAssigned() {
		int i;
		System.out.print("Assigned: [");
		for (i = 0; i < this.assigned.size(); i++)
			System.out.print(this.assigned.get(i)+ " ");
		System.out.println("]");
	}
	public void printRejected() {
		int i;
		System.out.print("Rejected: [");
		for (i = 0; i < this.rejected.size(); i++)
			System.out.print(this.rejected.get(i)+ " ");
		System.out.println("]");
	}
	
	public void printRejectNr() {
		int i;

		System.out.print("RejectNr: [");
		for (i = 0; i < this.rejectNumber.length; i++)
			System.out.print(this.rejectNumber[i]+ " ");
		System.out.println("]");
	}

	public int removeFromUnassigned() {
		
		return this.unassigned.remove(0);
	}

	public void replaceUnassigned(int[] unass) {
		this.unassigned.clear();
		for (int i = 0; i < unass.length; i++)
			this.unassigned.add(unass[i]);
	}
	
}
 