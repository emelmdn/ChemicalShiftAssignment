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
 * Public class to save the CSTable values and references.
 * The necessary operators for these reference values are implemented in this class.
 */

public class CSTable {
	int cSTableCount;
	private ArrayList<CSTableReference> references;
	
	
	public CSTable(){
		cSTableCount = 0;
		this.references = new ArrayList<CSTableReference>();
	}
	
	public CSTable(CSTable cs) {
		int i;
		this.references = new ArrayList<CSTableReference>();
		this.cSTableCount = cs.cSTableCount;
		for (i = 0; i <cs.references.size(); i++){
			CSTableReference ref = new CSTableReference(cs.references.get(i));
			this.references.add(ref);
		}	
	}

	public void setcSTableCount(int i){
		this.cSTableCount = i;
	}
	
	public int getcSTableCount(){
		return this.cSTableCount;
	}

	public void addReference(CSTableReference myRef) {
		this.references.add(myRef);
	}

	public void printCSTable(){
		int i;
		System.out.println("CSTABLE   "+ cSTableCount);
		for (i = 0; i < this.references.size(); i++){
			this.references.get(i).printReference();
		}
	}

	public boolean isNumeric(String input) { 
		try {  
			Double.parseDouble(input);  
		    return true;
		} 
		catch(NumberFormatException e){  
		    return false;  
		}  
	}
	
	/*
	 * This function returns reference of the specified residue's specified atom
	 * if type=1 then return reference, if type=2 return dev value
	 */
	
	public double getReference (String residueName, String atomName, int type){
		int i, j;
		boolean found = false;
		double reference = 0;
		double dev = 0;
		String atomNameQ;
		String lastChar;
		
		//1st check: check exactly the same name
		for (i = 0; i < this.references.size(); i++){ // check all references
			if (this.references.get(i).getResidueName().toString().equals(residueName.toString())){  // find the residue
				if (this.references.get(i).getAtomName().toString().equals(atomName.toString())){  // atom is found
					if (type == 1)
						reference =  this.references.get(i).getReference();
					else if (type == 2)
						dev =  this.references.get(i).getReferenceDev();
					found = true;
					break;
				}
			}
		}
		//2nd check: if the atom's same name does not exist, switch the first character with Q and check whether exist
		if (found == false && atomName.length()>1){
			atomNameQ = "Q" + atomName.substring(1);
			for (i = 0; i < this.references.size(); i++){ // check all references
					if (this.references.get(i).getResidueName().toString().equals(residueName.toString())){  // find the residue
						if (this.references.get(i).getAtomName().toString().equals(atomNameQ.toString())){  // atom is found
							if (type == 1)
								reference =  this.references.get(i).getReference();
							else if (type == 2)
								dev =  this.references.get(i).getReferenceDev();
							found = true;
							break;
						}
					}
				}
		}
		
		//3rd check: if last char is ", then replace it by '
		if (found == false && atomName.length()>1){
			lastChar = atomName.substring(atomName.length()-1,atomName.length());
			if (lastChar.equals("\"")){
				atomNameQ = "Q" +  atomName.substring(1, atomName.length()-1)+ "'";
				for (i = 0; i < this.references.size(); i++){ // check all references
					//	if (report) System.out.println("check the reference: ("+ this.references.get(i).getAtomName()+") in residue ("+this.references.get(i).getResidueName()+")");
						if (this.references.get(i).getResidueName().toString().equals(residueName.toString())){  // find the residue
							if (this.references.get(i).getAtomName().toString().equals(atomNameQ.toString())){  // atom is found
								if (type == 1)	
									reference =  this.references.get(i).getReference();
								else if (type == 2)
									dev =  this.references.get(i).getReferenceDev();
								found = true;
								break;
							}
						}
					}
			}
		}
		//4th check: remove the numbers one by one from the end
		if (found == false && atomName.length()>1){
			atomNameQ = "Q" + atomName.substring(1);

			int counter = 0;
			lastChar = atomNameQ.substring(atomName.length()-1,atomName.length());
			while (isNumeric(lastChar)){
				counter++;
				atomNameQ = "Q"+ atomNameQ.substring(1, atomName.length()-counter);
				for (i = 0; i < this.references.size(); i++){ // check all references
						if (this.references.get(i).getResidueName().toString().equals(residueName.toString())){  // find the residue
						if (this.references.get(i).getAtomName().toString().equals(atomNameQ.toString())){  // atom is found
							if (type == 1)
								reference =  this.references.get(i).getReference();
							else if (type == 2)
								dev =  this.references.get(i).getReferenceDev();
							found = true;
							break;
						}
					}
				}
				lastChar = atomNameQ.substring(atomName.length()-1-counter, atomName.length()-counter);
			}
		}

		if (type == 1)
			return reference;
		else if (type == 2)
			return dev;
		else
			return 0;
	}
	
}
