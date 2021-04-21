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
 * Public class to implement the dihedral angles of the NMR experiments.
 */
public class DihedralAngle {
	
	int number;
	String name;
	int[] atoms = new int[4]; // the number of 4 atoms that define the dihedral angle
	int lastAtom;

	public DihedralAngle(){

	}
	
	public DihedralAngle(int number, String name, int[] atoms, int lastAtom){
		int i;
		
		this.number = number;
		this.name= name;
		for (i = 0; i< 4; i++)
			this.atoms[i] = atoms[i];
		this.lastAtom = lastAtom;
	}
	
	public void printDihedralAngle(){
		int i;
		
		System.out.print(number+"\t"+name+"\t");
		for (i = 0; i < atoms.length; i++)
			System.out.print(atoms[i]+"\t");
		System.out.println(lastAtom);
	}
}
