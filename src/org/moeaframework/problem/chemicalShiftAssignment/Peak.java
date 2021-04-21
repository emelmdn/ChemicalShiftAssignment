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
 * Public class peak to save all measured peak and expected peak attributes and to apply all related functions. 
 */
public class Peak {
	
	private int dimension;
	private int peakFileNumber; // peak number in XEASY input file (for ex: HNCA_exp.peaks or HNCA.peaks file)
	private int index; //peak index in the expectedPeaks or measuredPeaks ArrayList of Spectrum.class, starting from 0
	private double[] chemicalShifts;
	private Atom[] atoms;
	private int colorCode;
	private String spectrumType;
	private double peakVolume; 
	private double peakVolumeError;
	private String integrationMethod;
	private int anInt;
	private boolean assigned = false;
	private int spectrumIndex;  //This peak's spectrum index on ChemicalShiftAssignment.allSpectra. The first one is zero.
	private String spectrumName; //Name of the spectrum that this peak belongs to.
	private int assignmentIndex = -1; //Index of this peak in the ChemicalShiftAssignment.assignment
	private double peakVolumeFromFile; 
	private boolean[] diagonal;
	private ArrayList<Integer> equalPeaks;
	private boolean equalsAllocated = false;
	private boolean isConnected;
	private int residueNumber;
	private int residueRange;
	

	public Peak(){
		this.assigned = false;
	}
	

	/*
	 * Constructor of the Peak class with the given input values. 
	 */
	public Peak(int dimension, int peakFileNumber, double[] chemicalShifts2, int colorCode, String spectrumType, double peakVolume, double peakVolumeError, String integrationMethod, int anInt, String[] atoms, String specName, int specIndex){
		int i,j, residueNumber;
		String atomName;
		
		this.dimension = dimension;
		this.assigned = false;
		this.peakFileNumber = peakFileNumber;
		this.chemicalShifts = new double[dimension];
		this.atoms = new Atom[dimension];
		this.equalPeaks = new ArrayList<Integer>();
		for (i = 0; i < dimension; i++){
			this.atoms[i] = new Atom();
		}
		for (i =0; i< dimension; i++){
			this.chemicalShifts[i] = chemicalShifts2[i];
		}
		this.colorCode = colorCode;
		this.spectrumType = spectrumType;
		this.peakVolume = peakVolume;
		this.peakVolumeFromFile = peakVolume; 
		this.peakVolumeError = peakVolumeError;
		this.integrationMethod = integrationMethod;		
		this.anInt = anInt;
		for (i =0; i< dimension; i++){
			atomName = atoms[i];
			this.atoms[i].setName(atoms[i]);
			residueNumber = Integer.parseInt(atomName.substring(atomName.indexOf(".")+1));
			this.atoms[i].setResidueNumber(residueNumber);
			for (j = 0; j < i; j++) {
				if (this.atoms[j].equals(atoms[i])) {
					this.diagonal[i] = true;
					this.diagonal[j] = true;
				}
			}
		}
		this.spectrumName = specName;
		this.spectrumIndex = specIndex;
		this.diagonal = new boolean[dimension];
		for (i = 0; i < dimension; i++)
			this.diagonal[i] = false;
		this.equalPeaks = new ArrayList<Integer>();
		equalPeaks.clear();
	}
	
	/*
	 * Copy constructor
	 */
	public Peak (Peak m){
		int i;	
		this.equalPeaks = new ArrayList<Integer>();
		this.dimension = m.dimension;
		this.assigned = m.assigned;
		this.peakFileNumber = m.peakFileNumber;
		this.index = m.index;
		this.chemicalShifts = new double[dimension];
		this.atoms = new Atom[dimension];
		for (i = 0; i < dimension; i++){
			atoms[i] = new Atom();
		}
		for (i = 0; i< dimension; i++){
			this.chemicalShifts[i] = m.chemicalShifts[i];
			this.atoms[i].setName(m.atoms[i].getName());
			this.atoms[i] = m.atoms[i];
		}
		this.colorCode = m.colorCode;
		this.spectrumType = m.spectrumType;
		this.peakVolume = m.peakVolume;
		this.peakVolumeFromFile = m.peakVolumeFromFile;
		this.peakVolumeError = m.peakVolumeError;
		this.integrationMethod = m.integrationMethod;
		this.spectrumIndex = m.spectrumIndex;
		this.spectrumName = m.spectrumName;
		this.assignmentIndex = m.assignmentIndex;
		this.diagonal = new boolean[dimension];
		for (i = 0 ; i < dimension; i++)
			this.diagonal[i] = m.diagonal[i];
		for (i = 0;  i< m.getEquivPeakSize(); i++)
			this.equalPeaks.add(m.getEquivPeaks().get(i));
	}
	
	/*
	 * This function returns the equivalent peaks of the current peak. 
	 */
	private ArrayList<Integer> getEquivPeaks() {
		return this.equalPeaks;
	}


	public void printPeak(){
		int i;
		System.out.println("========Print Peak========");
		System.out.println("peak number in file: " + peakFileNumber);
		System.out.println("dimension: " + dimension);
		System.out.println(" chemical shifts: ");
		for (i = 0; i < dimension; i++)
			System.out.print(chemicalShifts[i] + " ");
		System.out.println();
		System.out.print("atoms: ");
		for (i = 0; i < this.atoms.length; i++)
			System.out.print(this.atoms[i].getName()+" ");
		System.out.println();
		System.out.println("color code: " + colorCode);
		System.out.println("spectrum type: " + spectrumType);
		System.out.println("peak volume: "+ peakVolume);
		System.out.println("pSelect:" + peakVolumeFromFile);
		System.out.println("peak volume error: "+ peakVolumeError);
		System.out.println("integration method code: " + integrationMethod);
		System.out.println("anInt: "+anInt+"\t");
		for (i = 0; i < dimension; i++)
			System.out.print(atoms[i].getName() + " ");
		System.out.println();
		System.out.println ("spectrum name: "+ this.spectrumName);
		System.out.println("spectrum index: "+ this.spectrumIndex);
		System.out.println("assignment index: "+ this.assignmentIndex);
		System.out.print("Diagonal: ");
		for (i = 0; i < dimension; i++)
			System.out.print(this.diagonal[i] +" ");
		System.out.println();
		System.out.println("================");
	}
	
	/*
	 * This function prints the measured peak in the .peaks (XEASY) format line:     34   8.710  59.624 128.627 1 U   6.427E+00  0.000E+00 e 0     0     0     0
	 */
	public void printPeakLine(){
		int i;
		System.out.print("\t"+peakFileNumber+"\t");
		for (i = 0; i < this.dimension; i++)
			System.out.print(chemicalShifts[i] + "\t");
		System.out.print(colorCode +"\t"+spectrumType +"\t"+peakVolume+"\t"+peakVolumeError+"\t"+integrationMethod+"\t");
		System.out.print(anInt);
		for (i = 0; i < this.dimension; i++){
			System.out.print("\t"+ this.atoms[i].getName());
		}
		System.out.println();
	}

	public int getIndex(){
		return this.index;
	}
	
	public void setIndex (int i){
		this.index = i;
	}
	
	public Atom[] getAtoms(){
		return atoms;
	}

	public void setAtoms(Atom[] atoms){
		this.atoms = atoms;
	}
	
	/*
	 * This function returns the i. atom of the peak.
	 */
	public Atom getAtom(int i){
		if (i >= this.atoms.length)
			System.out.println("getAtom method of Peak tries to get index that is longer than its size!");
		return atoms[i];
	}
	
	/*
	 * This function returns the chemical shift values of the peak. 
	 */
	public double[] getChemicalShifts() {
		return chemicalShifts;
	}
	
	
	public void setIsConnected(boolean value) {
		this.isConnected = value;
	}
	
	public boolean getIsConnected() {
		return this.isConnected;
	}
	
	public void setResidueNumber(int value) {
		this.residueNumber = value;
	}
	
	public int getResidueNumber() {
		return this.residueNumber;
	}
	
	public void setResidueRange(int value){
		this.residueRange = value;
	}
	
	public int getResidueRange() {
		return this.residueRange;
	}
	
	public double getChemicalShift(int i){
		return chemicalShifts[i];
	}
	
	/*
	 * This function returns the chemical shift value of the givena atom.
	 */
	public double getChemicalShiftOfAtom(String atomName){
		int i;
		for (i = 0; i < this.dimension; i++){
			if (atomName == this.atoms[i].getName()){
				return this.chemicalShifts[i];
			}
		}
		return -1;
	}
	
	public void setChemicalShifts(double[] chemicalShifts) {
		this.chemicalShifts = chemicalShifts;
	}
	
	public int getDimension() {
		return dimension;
	}
	public void setDimension(int dimension) {
		this.dimension = dimension;
	}
	public String getSpectrumType() {
		return spectrumType;
	}
	public void setSpectrumType(String spectrumType) {
		this.spectrumType = spectrumType;
	}
	public int getColorCode() {
		return colorCode;
	}
	public void setColorCode(int colorCode) {
		this.colorCode = colorCode;
	}

	public double getPeakVolumeFromFile(){
		return this.peakVolumeFromFile;
	}
	
	public int getPeakFileNumber() {
		return peakFileNumber;
	}

	public void setPeakFileNumber(int peakFileNumber) {
		this.peakFileNumber = peakFileNumber;
	}

	public double getPeakProbability(){
		return peakVolume;
	}
	
	public double getPeakVolume() {
		return peakVolume;
	}

	public void setPeakVolume(double peakVolume) {
		this.peakVolume = peakVolume;
	}

	public double getPeakVolumeError() {
		return peakVolumeError;
	}

	public void setPeakVolumeError(double peakVolumeError) {
		this.peakVolumeError = peakVolumeError;
	}

	public String getIntegrationMethod() {
		return integrationMethod;
	}

	public void setIntegrationMethod(String integrationMethod) {
		this.integrationMethod = integrationMethod;
	}
	
	public boolean getAssigned(){
		return this.assigned;
	}
	
	public void setAssigned(boolean assigned){
		this.assigned = assigned;

	}
	
	public int getSpectrumIndex(){
		return this.spectrumIndex;
	}
	
	public void setSpectrumIndex(int specIndex){
		this.spectrumIndex = specIndex;

	}
	
	public String getSpectrumName(){
		return this.spectrumName;
	}
	
	public void setSpectrumName(String specName){
		this.spectrumName = specName;

	}
	
	public int getAssignmentIndex(){
		return this.assignmentIndex;
	}
	
	public void setAssignmentIndex(int assignmentIdx){
		this.assignmentIndex = assignmentIdx;

	}
	public void markAssigned(){
		this.assigned = true;
		for (int i = 0; i < dimension; i++)
			this.atoms[i].setAssigned(true);
	}
	
	public void markUnassigned(){
		this.assigned = false;
		for (int i = 0; i < dimension; i++)
			this.atoms[i].setAssigned(false);
	}



	public double getpExists() {
		return this.peakVolume;
	}



	public boolean[] getDiagonal() {
		return this.diagonal;
	}

	public void clearEqualPeaks() {
		this.equalPeaks.clear();
	}

	public int getEquivPeakSize() {
		return this.equalPeaks.size();
	}

	public boolean equalsAllocated() {
		return this.equalsAllocated;
	}

	public void addToEqualPeaks(int id) {
		if (isNewEqualPeak(id)) {
			this.equalPeaks.add(id);
		}
		this.equalsAllocated = true;
	}
	private boolean isNewEqualPeak(int peakId) {
		for (int i = 0; i < this.equalPeaks.size(); i++) {
			if (peakId == this.equalPeaks.get(i))
				return false;
		}
		return true;
	}

}
