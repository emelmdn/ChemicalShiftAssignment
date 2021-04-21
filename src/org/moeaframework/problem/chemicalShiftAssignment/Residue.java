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
 * The public class for the residue related attributes and functions. 
 */
public class Residue {
	
	String name;
	int angleCount;
	int atomCount;
	int firstAtomNumber;
	int lastAtomNumber;
	private ArrayList<DihedralAngle> dihedralAngles;
	private ArrayList<Atom> atoms;
	
	public Residue (){
		
		this.dihedralAngles = new ArrayList<DihedralAngle>();
		this.atoms = new ArrayList<Atom>();
		
	}
	
	/*
	 * The constructor of Residue class with the given input attributes. 
	 */
	public Residue(String name, int angleCount, int atomCount, int firstAtomNumber, int lastAtomNumber, ArrayList<DihedralAngle> dihedralAngles, ArrayList<Atom> atoms){
		
		this.dihedralAngles = new ArrayList<DihedralAngle>();
		this.atoms = new ArrayList<Atom>();
		int i;
		this.name = name;
		this.angleCount = angleCount;
		this.atomCount = atomCount;
		this.firstAtomNumber = firstAtomNumber;
		this.lastAtomNumber = lastAtomNumber;
		for (i = 0; i < dihedralAngles.size(); i++)
			this.dihedralAngles.add(dihedralAngles.get(i));
		for (i = 0; i < atoms.size(); i++)
			this.atoms.add(atoms.get(i));
	}

	/*
	 * Copy constructor
	 */
	public Residue(Residue re) {
		int i;
		
		this.dihedralAngles = new ArrayList<DihedralAngle>();
		this.atoms = new ArrayList<Atom>();
		
		this.name = re.name;
		this.angleCount = re.angleCount;
		this.atomCount = re.atomCount;
		this.firstAtomNumber = re.firstAtomNumber;
		this.lastAtomNumber = re.lastAtomNumber;
		for (i = 0; i < re.dihedralAngles.size(); i++)
			this.dihedralAngles.add(re.dihedralAngles.get(i));
		for (i = 0; i < re.atoms.size(); i++)
			this.atoms.add(re.atoms.get(i));
	}

	public void printResidue(){
		int i;
		System.out.println();
		System.out.println("RESIDUE\t"+name+"\t"+angleCount+"\t"+atomCount+"\t"+firstAtomNumber+"\t"+lastAtomNumber);
		for (i = 0; i < dihedralAngles.size(); i ++)
			this.dihedralAngles.get(i).printDihedralAngle();
		for (i = 0; i < atoms.size(); i++)
			this.atoms.get(i).printAtom();
			
	}
	
	public int getAtomCount(){
		return this.atomCount;
	}
	public Atom getAtom(int i){
		return atoms.get(i);
	}
	
	public String getName(){
		return this.name;
	}
}
