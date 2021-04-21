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

/*
 * Public class CSTableReference to save the CSTable values and operations
 */
public class CSTableReference {
	int idNumber;
	String residueName;
	String atomName;
	int notUsedNumber; 
	double referenceFrequency;
	double dev, third, fourth; 
	

	
	public CSTableReference(){
		
	}
	
	public CSTableReference(int id, String rName, String aName, int n, double ref, double sec, double thi, double fou){
		
		this.idNumber = id;
		this.residueName = rName;
		this.atomName = aName;
		this.notUsedNumber = n;
		this.referenceFrequency = ref;
		this.dev = sec;
		this.third = thi;
		this.fourth = fou;
	}
	
	public CSTableReference(CSTableReference ref){
		
		this.idNumber = ref.idNumber;
		this.residueName = ref.residueName;
		this.atomName = ref.atomName;
		this.notUsedNumber = ref.notUsedNumber;
		this.referenceFrequency = ref.referenceFrequency;
		this.dev = ref.dev;
		this.third = ref.third;
		this.fourth = ref.fourth;
	}

	public void printReference() {
		System.out.println(idNumber +"\t"+residueName+"\t"+atomName +"\t"+notUsedNumber+"\t"+referenceFrequency+"\t"+dev+"\t"+third+"\t"+fourth);
	}
	
	public double getReference(){
		return this.referenceFrequency;
	}

	public String getResidueName() {
		return this.residueName;
	}
	
	public String getAtomName(){
		return this.atomName;
	}

	public double getReferenceDev() {
		return this.dev;
	}
}
