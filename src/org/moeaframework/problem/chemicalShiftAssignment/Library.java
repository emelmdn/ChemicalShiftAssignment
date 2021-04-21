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
 * Public class Library to save and apply all  library attributes and operations. 
 */
public class Library {
	private int atomTypeCount;
	private ArrayList<AtomType> atomTypes;
	private ArrayList<Residue> residues;
	private CSTable cSTable;

	public Library(){
		this.atomTypes = new ArrayList<AtomType>();
		this.residues = new ArrayList<Residue>();
		cSTable = new CSTable();
	}
	public Library(int at){
		this.atomTypeCount = at;
		this.atomTypes = new ArrayList<AtomType>();
		this.residues = new ArrayList<Residue>();
		cSTable = new CSTable();
	}
	
	public Library (Library lib){
		int i;
		this.atomTypes = new ArrayList<AtomType>();
		this.residues = new ArrayList<Residue>();
		cSTable = new CSTable();
		
		this.atomTypeCount = lib.atomTypeCount;
		
		for (i = 0; i <lib.getAtomTypes().size(); i++){
			AtomType at = new AtomType(lib.getAtomType(i));
			this.atomTypes.add(at);
		}
		for (i = 0; i <lib.getResidueCount(); i++){
			Residue re = new Residue(lib.getResidue(i));
			this.residues.add(re);
		}	
		
		CSTable cs = new CSTable(lib.cSTable);
		this.cSTable = cs;
		
	}
	
	public void addAtomType(AtomType at){
		AtomType myAt = new AtomType (at);
		this.atomTypes.add(myAt);
	}
	
	public AtomType getAtomType(int i){
		return this.atomTypes.get(i);
	}
	
	public void addResidue(Residue re){
		Residue myRe = new Residue (re);
		this.residues.add(myRe);
	}
	
	public Residue getResidue(int i){
		return this.residues.get(i);
	}
	
	
	public ArrayList<AtomType> getAtomTypes(){
		return this.atomTypes;
	}
	
	public void setAtomtypes(ArrayList<AtomType> at){
		this.atomTypes = at;
	}
	
	public int getAtomTypeCount(){
		return atomTypeCount;
	}
	
	public void setAtomTypeCount(int at){
		this.atomTypeCount = at;
	}
	
	public int getResidueCount(){
		return this.residues.size();
	}
	public void printLibrary(){
		int i;
		System.out.println("**********************LIBRARY*********************");
		System.out.println("ATOMTYPES\t"+this.atomTypeCount);
		for (i = 0; i<this.atomTypes.size(); i++)
			this.atomTypes.get(i).printAtomTypeLine();
		for (i = 0; i < this.residues.size(); i++)
			this.residues.get(i).printResidue();
		this.cSTable.printCSTable();
		System.out.println("*******************************************");
	}
	
	public void setCSTableCount(int i){
		this.cSTable.setcSTableCount(i);
	}
	
	public int getcSTableCount(){
		return this.cSTable.getcSTableCount();
	}
	
	public void addReferenceToCSTable(CSTableReference ref){
		CSTableReference myRef = new CSTableReference(ref);
		this.cSTable.addReference(myRef);
	}
	
	/*
	 * This function assigns the reference frequency statistics from CSTABLE to each atom
	 */
	public void assignStatistics(ArrayList<String> aminoacidSequence){
		int i, j;
		String residueName;
		double referenceFrequency = 0;
		
		for (i = 0; i < this.residues.size(); i++){   //assign for each residue
			residueName = this.residues.get(i).getName();
			residueName = residueName.replaceAll("\\s+","");
			if (aminoacidSequence.contains(residueName)){ //residue exists in the sequence file
				for (j = 0; j < this.residues.get(i).getAtomCount(); j++){  // assign for each atom of that residue
					referenceFrequency = this.cSTable.getReference(this.residues.get(i).getName(), this.residues.get(i).getAtom(j).getName(),1);
					this.residues.get(i).getAtom(j).setReferenceFrequency(referenceFrequency);
				}
			}
		}
	}
	
	public double getReferenceFrequency (String residueName, String atomName){
	
		return this.cSTable.getReference(residueName, atomName, 1);
	}
	
	public double getReferenceDev (String residueName, String atomName){

		return this.cSTable.getReference(residueName, atomName, 2);
	}

}
