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
 * Public class ShiftReference to save the reference chemical shift values
 */
public class ShiftReference {

	private int atomNumber;
	float chemicalShift;
	float chemicalShiftError;
	String atomName;
	int residueNumber;
	
	public  ShiftReference() {
	
	}
	
	public  ShiftReference(int atomNumber, float chemicalShift, float chemicalShiftError, String atomName, int residueNumber) {
		this.atomNumber = atomNumber;
		this.chemicalShift = chemicalShift;
		this.chemicalShiftError = chemicalShiftError;
		this.atomName = atomName;
		this.residueNumber = residueNumber;
	}

	public ShiftReference(ShiftReference sr) {
		this.atomNumber = sr.atomNumber;
		this.chemicalShift = sr.chemicalShift;
		this.chemicalShiftError = sr.chemicalShiftError;
		this.atomName = sr.atomName;
		this.residueNumber = sr.residueNumber;
	}
	
	
}
