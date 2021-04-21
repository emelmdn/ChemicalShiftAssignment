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
 * Public class for the types of the atoms. 
 */
public class AtomType {
	
	int atomTypeNumber;
	String atomType;
	float coreRadius;
	int hydrogenCode; //1 for hydrogen atoms that can form hydrogen bonds, -1 for hydrogen bond acceptors (e.g. oxygens), and 0 for atoms that cannot be involved in hydrogen bonds
	int atomOrder; //0 for pseudo atoms, 1 for hydrogen, 6 for carbon, 7 for nitrogen 
	float lastValue; 
	
	public AtomType(){
		
	}

	public AtomType(int atomTypeNumber, String atomType, float coreRadius, int hydrogenCode, int order, float lastValue){
		this.atomTypeNumber = atomTypeNumber;
		this.atomType = atomType;
		this.coreRadius = coreRadius;
		this.hydrogenCode = hydrogenCode;
		this.atomOrder = order;
		this.lastValue = lastValue;
	}
	
	public AtomType(AtomType at) {
		this.atomTypeNumber = at.atomTypeNumber;
		this.atomType = at.atomType;
		this.coreRadius = at.coreRadius;
		this.hydrogenCode = at.hydrogenCode;
		this.atomOrder = at.atomOrder;
		this.lastValue = at.lastValue;

	}

	public void printAtomTypeLine(){
		System.out.println("\t"+atomTypeNumber+"\t"+atomType+"\t"+coreRadius+"\t"+hydrogenCode+"\t"+atomOrder+"\t"+lastValue);
	}
}
