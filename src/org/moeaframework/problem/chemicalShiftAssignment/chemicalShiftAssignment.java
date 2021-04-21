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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Random;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import javax.print.attribute.standard.PrinterMoreInfoManufacturer;

import org.moeaframework.core.NondominatedSortingPopulation;
import org.moeaframework.core.Problem;
import org.moeaframework.core.Solution;
import org.moeaframework.core.variable.BinaryIntegerVariable;
import org.moeaframework.core.variable.EncodingUtils;
import org.moeaframework.util.io.CommentedLineReader;

/*
 * The public class to implement new problem: ChemicalShiftAssignment
 * 
 */
public class chemicalShiftAssignment implements Problem {

	String path = "src/org/moeaframework/problem/chemicalShiftAssignment/";
	private ArrayList<Spectrum> allSpectra;
	private String[] spectrumNames;
	private String sequenceFile;
	private String shiftReferenceFile;
	private String libraryFile;
 	private ArrayList<String> aminoacidSequence;
 	private ArrayList<ShiftReference> shiftReferences;
 	private Library library;
 	private ArrayList<Atom> atoms;
 	private int[] expectedPeakSpecId;  //shows the specId of each expected peak of the whole list. (for ex: 28. exp peak belongs to HNCa, which is 0. spectrum
 	private int[] expectedPeakExpIdInSpec; //shows the expectedPeak Id in its own spectrum of each expected peak of the whole list. (For ex: 28. exp peak is 3. exp Peak in its own spectrum HNCA)
	
 	private int nObjectives; //is set in constructor below
 	private int nVariables; // is set to expectedPeak number of first spectrum for MOEA in constructor below
 	private int nConstraints; // is set to 1, but only the first one is used; others are constant as of now
 	private double weightTwo = 2;
 	private double maximumScore;
	private double maximumDeviation;	
	private static double smallestValue = 1.0E-30;
	private double weightThree = 0.49682104549742656;
	private double currentG = 0; //used in evaluateObjectives method to save the G value of the current assignment
 	private double factor = 3.0; 
 	long seed = 428; 
	Random r = new Random (seed);
	private double gaussian; 
	private boolean gaussianStored  = false;
	private double NNN = 999;
	private Assignment assignment;
	private int residueAssignment = 6;
	private int  residueExpression = 8;   
	private double[] assignmentDistance = {1.0,0.99,0.98,0.97,0.96,0.5,0.001};
	private double[] expScore = {1.0,0.98,0.95,0.90,0.8,0.6,0.5,0.15,0.001};
	private double factorThreshold = 1;
	private double unassignedRatio = 0.2;
 	private int assignmentType; 
 	private double averageAtom = 0;
 	private double tempValue = 0.2; 
 	private int addFrequencyCounter = 0;
 	private int numberOfSpectrumTypes = 1; 
 	private double[] spectrumWeight = new double[1];
 	private int weightTypeValue = 0; 
 	private boolean[] holdAssignment; 
 	private double degeneracyFactor = 1.0;
 	private int iopt = -1; 
 	private int thisWeightOne = 4; 
 	private double thisMaximumScore = 0;
 	
 	private double[] populationScore;
 	private int[] populationActive;
 	private double temp = 0;
 	NondominatedSortingPopulation population;
 	double probabilityDegeneracy;
 	private Sequence sequence;
 	private int populationSize;
 	private double[] maximum;
 	private double[] minimum;
 	private double[][] middleScore;
 	private double[][] scoreResidue;
 	private double[][] populationConnectAtom;
 
 	
 	/*
 	 * The constructor that ready the input files and applies the initialization to the whole problem model attributes accordingly.
 	 * 
 	 */
	public chemicalShiftAssignment(File file) throws IOException {
		super();
		this.allSpectra = new ArrayList<Spectrum>();
	
		this.aminoacidSequence = new ArrayList<String>();
		this.shiftReferences = new ArrayList<ShiftReference>();
		this.atoms = new ArrayList<Atom>();
		int assignmentSize, i;
		this.library = new Library();
		this.population = new NondominatedSortingPopulation();
		boolean report = true;
		
		loadFileNames(new FileReader(file));
		this.printConfigurationFile();
		loadFiles();
		if (report) System.out.println("All files are loaded successfuly.");
		initialize();
		if (report) System.out.println("Initialization is done successfuly.");
		assignmentSize = 0;
		for (i = 0; i< this.allSpectra.size(); i++)
			assignmentSize += this.allSpectra.get(i).getExpectedPeaksSize();
		if (report) System.out.println("Assignment size: "+ assignmentSize);
		this.assignment = new Assignment(assignmentSize);
	
		this.nVariables = this.getExpectedPeakSize();
		setObjectiveConstraintCounts();  //constraint and objective counts are updated 
	
	}
	
	/*
	 * This function iterates the hill climber algorithm. 
	 */
	public Solution runHillClimber(Solution parent) {
		int j;
		Solution solution = this.newSolution();
		this.resetAssignment();
		this.initializePeaks();
	
		fillAssignmentWithCurrentSolution(parent);
		
		iterateHillClimber();
			
		for (j = 0; j < solution.getNumberOfVariables(); j++) {
			BinaryIntegerVariable variable = (BinaryIntegerVariable)solution.getVariable(j);
			variable.setValue(this.assignment.getAssignment(j));
		}
		
		return solution;
		
	}
	
	/*
	 * This function implements the iteration of the hill climber algorithm. 
	 */
	private void iterateHillClimber () {
		
		int i, counter = 10000;
		
		
		for (i = 0;  i < counter; i++) {
			updateAssignment();
			makeAssignment();
			
		}
		
	}
	
	/*
	 * This function updates the assignment with the current solution 
	 */
	private void fillAssignmentWithCurrentSolution(Solution parent) {
		int expId, measId;
		
		for (expId = 0; expId < parent.getNumberOfVariables(); expId++) {
		
			
			BinaryIntegerVariable variable = (BinaryIntegerVariable)parent.getVariable(expId);
			measId = variable.getValue();
			
			if (measId > -1) { //expected peak is assigned to a measured peak
				 assignExpPeakToMeasPeak(expId, measId);
			} else if (measId == -1) {
				rejectPeakAssignment(expId);
			}
		}
		
	}

	/*
	 * This function creates one individual from scratch and sends the assignment as Solution to be used in the initialization of population
	 */
	protected Solution initializeOneIndividual(int i) {
		
		this.iopt = i+1; 
		
		Solution solution = this.newSolution();
		
		resetAssignment(); 

		solution = initializeIndividual();
		return solution;
	}
	
	/*
	 * This function initializes the list with all expected peaks as unassigned and randomized.
	 */
	protected void initializePeaks() {
		int unassignedSize, i;
		int[] positions;
		int[] values;
		
		unassignedSize = this.assignment.getUnassigned().size();
		if (unassignedSize> 0) {
			positions = new int[unassignedSize];
			values = new int[unassignedSize];
			
			for (i = 0; i < unassignedSize; i++) {
				positions[i] = i;
				values[i] = this.assignment.getUnassigned().get(i);
			}
			
			this.quickSort(values);
			this.randomizeArrayValues(values);
			this.assignment.replaceUnassignedWith(values);
			
		}
	}
	
	/*
	 * This function initializes one individual and returns it.
	 */
	public Solution initializeIndividual(){
	
		int expectedPeak = -1;
		Solution solution = this.newSolution();
		this.assignmentType = 1; 
		
	
		initializePeaks();
		makeAssignment();

		//first assignment of the child 
		this.assignment.resetUnassigned();
		while(true){ //find rejected peaks and assign them to unassigned list
			expectedPeak = getRejectedExpectedPeak();
			if (expectedPeak == -1)
				break;
			this.assignment.addToUnassigned(expectedPeak);
		}
		makeAssignment();
		assignPeaks();
		
		calculateScore(); 
		for (int j = 0; j < solution.getNumberOfVariables(); j++) {
			BinaryIntegerVariable variable = (BinaryIntegerVariable)solution.getVariable(j);
			variable.setValue(this.assignment.getAssignment(j));
		}
		return solution;
		
	}
	
	/*
	 * This function updates the score values according to the current assignment.
	 */
	public void saveAssignment() {
		int i;
		int position = this.iopt -1;
		
		this.populationScore[this.iopt-1] = this.currentG;
		for (i = 0; i < this.sequence.getNumber(); i++) {
			this.middleScore[position][i] = this.sequence.getScoreResidue(i);
			this.scoreResidue[position][i] = this.sequence.getScore(i);
			this.populationConnectAtom[position][i] = this.sequence.getConnectSequence(i);
		}
	}
	
	/*
	 * This function calculates the score of the current assignment. 
	 */
	private void calculateScore() {
		int atomId, i;
		double preValue;		
		int popSize = this.population.size();
		
	
		for (i = 0; i < popSize; i++)
			this.populationScore[i] = 0;
		
		for (atomId = 0 ; atomId < this.atoms.size(); atomId++) {
			if (this.atoms.get(atomId).getExpectedPeakSize()> 0) { 
			
				if (this.atoms.get(atomId).getPhaseNr() < 10) {
					preValue = 0.8 + this.atoms.get(atomId).getPhaseNr()*0.02;
				} else {
					preValue = 0.99;
				}
				this.atoms.get(atomId).setAverageScoreLocal(preValue * this.atoms.get(atomId).getAverageScoreLocal() + (1.0-preValue) * this.atoms.get(atomId).getScoreLocal());
		        this.atoms.get(atomId).setDeviationLocal(preValue * this.atoms.get(atomId).getDeviationLocal() + (1.0-preValue) * Math.pow(this.atoms.get(atomId).getAverageScoreLocal() - this.atoms.get(atomId).getScoreLocal(), 2));
				
			}
			        			
		}
	}
	
	/*
	 * This function triggers the constructive initialization method for the initial population creation.
	 */
	private void assignPeaks(){
		int runs = 30000;
		int irun;
		
		
		for (irun = 1; irun <= runs; irun++){
				
			updateAssignment();
			makeAssignment();
		}
	}
	
	/*
	 * This function searches for the expected peaks to be assigned for a better solution during population initialization.
	 */
	private void updateAssignment() {
		double rand;
		int expPeakId = -1, spectrumId;
		int peakId = -1,dimension;
		boolean [] changeAtom;
	
		rand = this.getRandomNumber();
		
		rand = rand * (this.unassignedRatio + 1/this.unassignedRatio);
		if (rand < this.unassignedRatio){ 
			expPeakId = getDegeneratePeak();
		}
		if (expPeakId == -1){
			expPeakId = getRejectedExpectedPeak();
		}
		if (expPeakId >= 0){
			
			if (this.assignmentType == 2) {  //evolOpt
	            peakId = getEvolutionaryPeakId (expPeakId);
		         
	                
			}  else {  //localOpt
	            peakId = getLocalPeakId (expPeakId);
			}
			 if (peakId != -1) { 
	        	spectrumId = this.expectedPeakSpecId[expPeakId];
	    		dimension = this.allSpectra.get(spectrumId).getDimension();
	    		changeAtom = new boolean[dimension];
	    		changeAtom = checkAtom (expPeakId, peakId);
	            resetPeaks (expPeakId, changeAtom);
	            
	            assignExpPeakToMeasPeak(expPeakId, peakId);
	 
	        } else{
	        	if (this.assignment.getStatus(expPeakId) == "PR"){ //if rejected       		
	        		
	        		this.assignment.addToUnassigned(expPeakId);
	        		this.assignment.setAssignment(expPeakId, -1);
	        		
	        	}else if (this.assignment.getStatus(expPeakId) == "PA"){ //if accepted
	    	        
	        		
	        		this.assignment.addToAssigned(expPeakId);
	        	}
	        }
		}
		
     
	}
	
	/*
	 * This function sends the expected peak indices (global expId, so not spec specific) that are mapped to the measId in the same spectrum
	 */
	private ArrayList<Integer> getAssignmentPeaks(int measuredPeakId, int spectrumId) {
		int i, startIndex, endIndex;
		ArrayList<Integer> assignmentPeaks = new ArrayList<Integer>();
		
		if (this.allSpectra.get(spectrumId).getExpectedPeaksSize() > 0 ) {
			startIndex = this.allSpectra.get(spectrumId).getAssignmentIndexFirstExpPeak();
			endIndex = this.allSpectra.get(spectrumId).getAssignmentIndexLastExpPeak();
			for (i = startIndex; i <= endIndex; i++) {
				if (this.assignment.getAssignment(i) == measuredPeakId) {
					assignmentPeaks.add(i);
				}
			}
		}
		return assignmentPeaks; 
	}
	
	/*
	 * This function sends the same values as in the arraylist taken as input but in a different variable
	 */
	private ArrayList<Integer>  getList(ArrayList<Integer> assignList) {
		ArrayList<Integer> linkList = new ArrayList<Integer>();
		int i;
	
		for (i = 0; i < assignList.size(); i++)
			linkList.add(assignList.get(i));
		
		return linkList;
		
	}
	/*
	 * This function returns expected peak which is part of ambiguous assignment
	 */
	private int getDegeneratePeak() {
		double randomNumber, threshold;
		int pExpPeak = -1, expPId = -1;
		int i, assignmetn, spectrumId, expIdInSpec, assignedPeakCounter, degenerateAssignmentCounter;
		ArrayList<Integer> assignmentArrayList = new ArrayList<Integer>(); 
		ArrayList<Integer> linkList = new ArrayList<Integer>();
				
		randomNumber = this.getRandomNumber();
	
		threshold = this.degeneracyFactor * randomNumber;
		
		i = 0;
		
		proceed: while (true) { 
			//find first expPeak in assignmentList with status assigned
			while (true) {
				pExpPeak = -1;
				if (i == 100)
					break proceed;
				expPId = this.assignment.dequeueFromAssignedList();
				
				if (expPId == -1) {
					break proceed;
				}
				pExpPeak = expPId;
				if ( this.assignment.getStatus(pExpPeak) == "PA") {
					break;
				}
			}
			assignmetn = this.assignment.getAssignment(pExpPeak);
			spectrumId = this.expectedPeakSpecId[pExpPeak];
			expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
			assignedPeakCounter = this.allSpectra.get(spectrumId).getAssignedPeakCount(assignmetn);
			if (assignedPeakCounter == 1) {
				degenerateAssignmentCounter = 0;
			} else {
				assignmentArrayList = getAssignmentPeaks(assignmetn, spectrumId);  
				linkList = getList(assignmentArrayList);
				degenerateAssignmentCounter = getNumberOfDegenerateAssignments(linkList,expPId);
				
			}
			if (assignedPeakCounter != (degenerateAssignmentCounter+1) && (this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).getPeakVolumeFromFile() >= threshold)) {
				break proceed;
			}
			this.assignment.addToAssigned(expPId);
			threshold = threshold * 0.8;
			i += 1;
			
		}
		
		if (i > 10) {
			this.degeneracyFactor = this.degeneracyFactor * 0.8;
		} else if (i < 3 && expPId != -1) {
			this.degeneracyFactor = this.degeneracyFactor * 1.2;
		}
		return pExpPeak;
	}
	
	/*
	 * This function calculates the total number of degenerate assignments and returns this value as output.
	 */
	private int getNumberOfDegenerateAssignments(ArrayList<Integer> peaks, int peakId) {
		int degenerateCounter, spectrumId, expIdInSpec, expIdInSpecAssId, d, iterator, assignmentId, dim, atomId, i;
	
	
		spectrumId = this.expectedPeakSpecId[peakId];
		expIdInSpec = this.expectedPeakExpIdInSpec[peakId];
		Peak expectedPeak = new Peak();
		ArrayList<Atom> atoms = new ArrayList<Atom>();
		
		d = this.allSpectra.get(spectrumId).getDimension();
	

		
		for (i = 0; i < d; i++)
			atoms.add(this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).getAtom(i));
	
		degenerateCounter = 0;
		
		for (iterator = 0; iterator < peaks.size(); iterator ++) {
	
			assignmentId = peaks.get(iterator);
			expIdInSpecAssId = this.expectedPeakExpIdInSpec[assignmentId];
			spectrumId = this.expectedPeakSpecId[assignmentId]; 
			d = this.allSpectra.get(spectrumId).getDimension();
			if (peakId != assignmentId) {
				
				expectedPeak = this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpecAssId);
				dim = 0;
				for (atomId = 0; atomId < d; atomId++) {
					if ((expectedPeak.getAtom(atomId).getName().equals(atoms.get(atomId).getName()))) {
						dim += 1;
					}
				}
				if (dim == (d-1)) {
					degenerateCounter += 1;
				}
			}
		}
		
		return degenerateCounter;

	}

	/*
	 * This function resets the whole assignment and sets all expected peaks to unassigned status.
	 */
	private void resetPeaks(int expPeaks, boolean[] save) {
		
		int spectrumDimension, spectrumId, atomId, atom, expIdInSpec, atomsIndexInChemShift, atomExpectedId, peakId;

		spectrumId = this.expectedPeakSpecId[expPeaks];
		spectrumDimension = this.allSpectra.get(spectrumId).getDimension();
		
		resetPeakAssignment(expPeaks);
		atomId = -1;
	
	    dimension: while(true){
	    	while(true){   
	    		atomId = atomId +1;
	    		if (atomId >= spectrumDimension){
	    			break dimension;
	    		}
	    		if (!(save[atomId])){
	    			break;
	    		}
	    	}
	    	
	    	expIdInSpec = this.expectedPeakExpIdInSpec[expPeaks];
			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(atomId).getAtomsIndex();
			
			atom = atomsIndexInChemShift;
			startPeakId(atom);
			
	    	
			while (true){ 
				peakId = getNextPeak(atom, false);
				if (peakId == -1){
					break;
				}
			    atomExpectedId = peakId;
			  
	            if (this.assignment.getStatus(atomExpectedId) != "PU"){ 
	            
	            	resetPeakAssignment(atomExpectedId);
	            	
	            	this.assignment.addToUnassigned(peakId); //add çaðýrmam lazýmdý, methodu yazdým ama burada da böylece hallettim
	            	this.assignment.setExpStatusToUnassigned(peakId);
	            	
	            } 
			}
	    }


	}

	/*
	 * This function deletes the assigned expected peaks of the given atom.
	 */
	private void startPeakId(int atomId){
		this.atoms.get(atomId).nullifyiEP();
	}

	/*
	 * This function is called from resetPeakASsignments to reset all existing peak assignments. 
	 */
	private void resetPeakAssignment(int pExpPeak){
		int spectrumId, dimension, assignmentId, specId, dim, expIdInSpec, atomsIndexInChemShift, atomIndexInSpec, atomId;
		double chemShift;
		
	
		spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();


		if (this.assignment.getAssignment(pExpPeak)> -1){
	      
			assignmentId = this.assignment.getAssignment(pExpPeak);
			specId = this.expectedPeakSpecId[pExpPeak];
			
			for (dim = 0 ; dim < dimension; dim ++){
				chemShift = this.allSpectra.get(specId).getChemicalShiftValues(dim, assignmentId);
				expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
				atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(dim).getAtomsIndex();
				atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(dim).getIndex();
				atomId = atomsIndexInChemShift;
		
				 deleteFrequency(atomId, chemShift, specId, dim, atomIndexInSpec);
		            
			}
			removeConnection(pExpPeak);
	        
		}
		
		this.assignment.setAssignment(pExpPeak, -1);
	    
	    	    
	}
	
	/*
	 * This function is called from the resetPeakAssignment and it deletes the given expected peaks assignments from all related objects, such as assignment and spectrum.
	 */
    private void removeConnection(int pExpPeak) {
    	int spectrumId, dimension, assignmentId, assignmentCount;
    
    	spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		assignmentId = this.assignment.getAssignment(pExpPeak);
		this.assignment.setAssignment(pExpPeak, -1);
		assignmentCount = this.allSpectra.get(spectrumId).getAssignedPeakCount(assignmentId);
    	this.allSpectra.get(spectrumId).setAssignedPeakCount(assignmentId, assignmentCount -1 );      
    
	}

    private void setAverageAtom(double average){
    	this.averageAtom = average;
    }

    /*
     * This function removes the given frequency from the atom. 
     */
	private void deleteFrequency(int atomId, double frequency, int spectrumId, int idim, int atomIndexInSpec){
    	int numerator;
    	double total;
		
		this.atoms.get(atomId).setAssignedMeasuredPeakCount(this.atoms.get(atomId).getAssignedMeasuredPeakCount() -1);
		if (this.atoms.get(atomId).getAssignedMeasuredPeakCount() == 0){
			resetMeasuredFrequency(atomId, spectrumId, atomIndexInSpec);
		}else {  
			this.atoms.get(atomId).setMeasPFreqSum(this.atoms.get(atomId).getMeasPeakFreqSum() - frequency);
			total = this.atoms.get(atomId).getMeasPeakFreqSum();
			numerator = this.atoms.get(atomId).getAssignedMeasuredPeakCount();
			this.atoms.get(atomId).setMeasuredPeakAverage(total / numerator);
			numerator = this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).getAssignedPeakCountSpectrum();
			this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setAssignedPeakCountSpectrum(numerator-1);
			if (this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).getAssignedPeakCountSpectrum() == 0){
				this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setMeasuredFrSpectrum(0);
				this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setMeasuredPeakFreqSum(0);
				
			}else{
				total = this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).getMeasuredPeakFreqSumSpectrum();
				this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setMeasuredPeakFreqSum(total -frequency);
				total = this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).getMeasuredPeakFreqSumSpectrum();
				numerator = this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).getAssignedPeakCountSpectrum();
						
				this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setMeasuredFrSpectrum(total / numerator);
			}
		}
		this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setScored(false); 
		this.atoms.get(atomId).setScored(false);
	

	}
	
	
	/*
	 * This function is used to reset the measured chemical shift frequency of the given atom in the given spectrum
	 */
	private void resetMeasuredFrequency(int atomId, int spectrumId, int atomIndexInSpec) {
	
	    
		this.atoms.get(atomId).resetMeasuredFrequencyValu();
	    
		this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).resetMeasFSpec();
		this.atoms.get(atomId).setAssigned(false); 
		this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setAssigned(false);;
		this.atoms.get(atomId).increasecf2ByOne();
		this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).increasecf2ByOne();
		this.atoms.get(atomId).setScored(false);
		this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setScored(false);
	}

	/*
	 * This function checks the atoms' chemical shift values for the given expected peak and
	 * sends true if they should be improved
	 * sends false if they are already in the predefined deviation boundaries.
	 */
	private boolean[] checkAtom(int pExpPeak, int peakId){
		int spectrumId, dimension, iterator, expIdInSpec, atomsIndexInChemShift, atomIndexInSpec, atomId;
		double measuredFrequency, deviationMeasuredFrequency;
		boolean [] adapt;
		
		
		spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		adapt = new boolean[dimension];
      
		for (iterator = 0 ; iterator < dimension; iterator ++){
			
			expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			atomId = atomsIndexInChemShift;
			if (this.atoms.get(atomId).getAssigned(this.assignment.getAssignmentWhole()) == true){
				getMeasuredFrequency(spectrumId, expIdInSpec, iterator); 
				measuredFrequency = this.atoms.get(atomId).getMeasFreqDynamic();
				
				deviationMeasuredFrequency = this.atoms.get(atomId).getMeasDevDynamic(); 
				if (this.allSpectra.get(spectrumId).getChemicalShiftValues(iterator, peakId) >= (measuredFrequency - deviationMeasuredFrequency) &&
						this.allSpectra.get(spectrumId).getChemicalShiftValues(iterator, peakId) <= (measuredFrequency + deviationMeasuredFrequency)){
					adapt[iterator] = true;
				}else{
					adapt[iterator] = false;
				}
			} else {
				adapt[iterator] = true;
			}
		}
        
    	return adapt;
	}

	
	/*
	 * This function searches for a measured peak that should be assigned to the given expected peak during the population initialization constructive logic. 
	 */
	private int getEvolutionaryPeakId(int pExpPeak) {
		int peakId = -1, spectrumId, dimension, iterator, atomId, atomsIndexInChemShift, atomIndexInSpec, expIdInSpec;
		boolean [] assignment;
		double random;
		
		spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		this.holdAssignment = new boolean[dimension];
		
		assignment = new boolean[dimension];
	    if (replaceAtom(pExpPeak) > 0){   
	
	    	if (true){ 
	    		spectrumId = this.expectedPeakSpecId[pExpPeak];
	    		for (iterator =0 ; iterator<dimension; iterator++){
	    			expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
	    			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
	    			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
	    			atomId = atomsIndexInChemShift;
	    			assignment[iterator] = this.atoms.get(atomId).getAssigned(this.assignment.getAssignmentWhole());
	    			getFrequencyRange (spectrumId, expIdInSpec, iterator, true, holdAssignment[iterator]);
	    			getExpectedFrequency(spectrumId, expIdInSpec, iterator);
	    			getMeasuredFrequency(spectrumId, expIdInSpec, iterator);
	            }
	    		expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
	    		peakId= selectedEvolutionary(spectrumId, expIdInSpec, holdAssignment,holdAssignment);
	    		
	    		if (peakId == -1 && this.temp > 0) {
					random = getRandomNumber();
					random = Math.pow(random, this.temp);
					if (random <= 0.5) {
						if (this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).equalsAllocated()) {
							peakId = getEqualPeakId(spectrumId, expIdInSpec, assignment);
						}

						boolean assign = false;
						boolean diagonal = false;
						for (int k = 0; k < dimension; k++) {
							assign = assign || assignment[k];
							diagonal = diagonal ||  this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).getDiagonal()[k];
						}
						
		                if (peakId == -1 && (! (!assign && diagonal)) && random <= 0.3) {
		                	peakId = selectPeak(spectrumId, expIdInSpec, assignment);
		                }
					}
				}
        

	    		 
	    	}
	    }
	   	  
	    return peakId;
	}

	/*
	 * This function searches for a measured peak that should be assigned to the given expected peak during the improvement process of the initialization of population.
	 */
	private int getLocalPeakId(int pExpPeak) {
		int peakId = -1, spectrumId, dimension, iteraotr, atomId, atomsIndexInChemShift, atomIndexInSpec, expIdInSpec;
	
		boolean [] assignment;
		
		spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		this.holdAssignment = new boolean[dimension];
		
		assignment = new boolean[dimension];
	    if (replaceAtom(pExpPeak) > 0){  
	    	if (true){ 
	    		spectrumId = this.expectedPeakSpecId[pExpPeak];
	    		for (iteraotr =0 ; iteraotr<dimension; iteraotr++){
	    			expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
	    			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iteraotr).getAtomsIndex();
	    			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iteraotr).getIndex();
	    			atomId = atomsIndexInChemShift;
	    			assignment[iteraotr] = this.atoms.get(atomId).getAssigned(this.assignment.getAssignmentWhole());
	    			getFrequencyRange (spectrumId, expIdInSpec, iteraotr, true, holdAssignment[iteraotr]);
	    			getExpectedFrequency(spectrumId, expIdInSpec, iteraotr);
	    			getMeasuredFrequency(spectrumId, expIdInSpec, iteraotr);
	            }
	    		expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
	    		peakId= selectPeak(spectrumId, expIdInSpec, holdAssignment);
	    		 
	    	}
	    }
	    	
	  
	    return peakId;
	}

	
	/*
	 * This function checks the atoms' chemical shift values of the given expected peak and 
	 * returns true if there is improvement needed
	 * returns false otherwise.
	 */
	private int replaceAtom(int pExpPeak) {
		int spectrumId, dimension, switchLocal, iterator, atomId, expIdInSpec,atomsIndexInChemShift, atomIndexInSpec, switchAtom;
		double buffer = 0, average, offset, value;
		double [] score, meanValue, deviation;

		
		spectrumId = this.expectedPeakSpecId[pExpPeak];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		expIdInSpec = this.expectedPeakExpIdInSpec[pExpPeak];
		switchLocal = 0; 
		average = this.averageAtom;
		score = new double[dimension];
		meanValue = new double[dimension];
		deviation = new double[dimension];
		
		for (iterator = 0; iterator < dimension; iterator ++){
			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			atomId = atomsIndexInChemShift;
        	atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			atomId = atomsIndexInChemShift;
			updateScoreLocal(atomId);
			score [iterator] = this.atoms.get(atomId).getScoreLocal();
			meanValue [iterator] = this.atoms.get(atomId).getAverageScoreLocal();
			deviation [iterator] = this.atoms.get(atomId).getDeviationLocal();
			buffer = buffer + this.atoms.get(atomId).getScoreLocal();
		} 
		
        for (iterator = 0 ; iterator < dimension; iterator++){
        	atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			atomId = atomsIndexInChemShift;
        	score[iterator] = score[iterator] - meanValue[iterator];  
	        score[iterator] = score[iterator] / Math.max(Math.sqrt(deviation[iterator]), chemicalShiftAssignment.smallestValue);
			
        }
        for(iterator = 0 ; iterator < dimension; iterator++){
 
            if (average < score[iterator]){
                offset = 0.05 * this.tempValue; //0.01
            } else {   
                offset = 0.05 - 0.05 * this.tempValue;  //0.04
            }
            average = (1.0 - offset) * average + offset * score[iterator];
			
        }
        this.setAverageAtom(average);
     
        for(iterator = 0 ; iterator < dimension; iterator++){
            value = score[iterator] - average;
            atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			atomId = atomsIndexInChemShift;
			if (this.atoms.get(atomId).getAssigned(this.assignment.getAssignmentWhole()) == true) {
				holdAssignment[iterator] = true;
			} else { 
				holdAssignment[iterator] = false;
			}
		
			if (this.atoms.get(atomId).getAssigned(this.assignment.getAssignmentWhole()) == true && value < 0 ) { 
				holdAssignment[iterator] = false;
				switchLocal += 1;
			}
        }
        switchAtom = switchLocal;
	 	return switchAtom;
	}


	/*
	 * This function is called from the replaceAtom function to update the existing score values of the given atom.
	 */
	private void updateScoreLocal(int atomNumber){
	
		int peakId = -1, diagonal, assignmentId, spectrumId, peakIdInSpec, dimension, iterator, phaseNr;
		boolean found;
        double peakExists, vol, dist, DSIGMA = 0.8, distfaktor, oldValue;
        double scoreFactor;
        
        this.atoms.get(atomNumber).setScored(true);
        this.atoms.get(atomNumber).setLocalScore(0);
       
        if (this.atoms.get(atomNumber).getMeasuredPeakAverage() != 0){ 
        	peakId = getNextPeak(atomNumber,true);
        	while(peakId != -1) {
        		diagonal = 1;
        		assignmentId = this.assignment.getAssignment(peakId);
                if (assignmentId != -1){
                   spectrumId = this.expectedPeakSpecId[peakId];
                   peakIdInSpec = this.expectedPeakExpIdInSpec[peakId];
                   peakExists = this.allSpectra.get(spectrumId).getExpectedPeaks().get(peakIdInSpec).getpExists();
                   dimension = this.allSpectra.get(spectrumId).getDimension();
                   found = false;
                   for (iterator = 0; iterator < dimension; iterator++){
                	   if (this.allSpectra.get(spectrumId).getExpectedPeak(peakIdInSpec).getAtom(iterator).getAtomsIndex() == atomNumber){
                		   found = true;
                		   break;
                	   }
                   }
                   
                   if (found) {
                	
                	   peakExists = peakExists / this.allSpectra.get(spectrumId).getAssignedPeakCount(assignmentId);
                	   vol = 4.1809802569675831;  
                	   dist = this.allSpectra.get(spectrumId).getDValue(); 
                	 
                	   if (vol > 0 &&  dist > vol){
                		   dist = this.allSpectra.get(spectrumId).getDValue(); 
                           distfaktor = Math.exp(-0.5 * Math.pow((dist-vol)/DSIGMA, 2));
                           peakExists = peakExists * distfaktor;
                       }
                     
                	   scoreFactor = this.atoms.get(atomNumber).getScalingFactorSpectrum( this.allSpectra.get(spectrumId).getCategory()); 
                	   this.atoms.get(atomNumber).setLocalScore(this.atoms.get(atomNumber).getScoreLocal() + peakExists * scoreFactor);
                       
                   }
                }
                peakId = getNextPeak(atomNumber, false);
        	}
        }
        phaseNr = this.atoms.get(atomNumber).getPhaseNr();
        if (phaseNr < 10){
        	this.atoms.get(atomNumber).setPhaseNr(phaseNr +1);
        	oldValue =0.8 + this.atoms.get(atomNumber).getPhaseNr()*0.02;
        	this.atoms.get(atomNumber).setAverageScoreLocal(oldValue *this.atoms.get(atomNumber).getAverageScoreLocal() + (1 - oldValue) * this.atoms.get(atomNumber).getScoreLocal() ); 
        	this.atoms.get(atomNumber).setDeviationLocal(oldValue * this.atoms.get(atomNumber).getDeviationLocal() +(1 - oldValue) * Math.pow(this.atoms.get(atomNumber).getAverageScoreLocal()-this.atoms.get(atomNumber).getScoreLocal(), 2));
       }
	}
	
	
	/*
	 * This function sends the next expected peak of the given atom Id. 
	 */
	private int getNextPeak(int atomId, boolean start) {
		int id = -1;
	    
	    
	    if (this.atoms.get(atomId).getExpectedPeakSize()>0){
	    	
	            if (!(this.atoms.get(atomId).getiEPSize()>0) || start) { 
	                this.atoms.get(atomId).initializeiEP(this.atoms.get(atomId).getExpectedPeaks());
	                id = getId (atomId);
	            } else {
	                id = next (atomId);

	            }
	    }else {
	          id = -1;
	    }

		return id;
	}

	/*
	 * This function updates the expected peak list of the given atom and returns the next expected peak Id of the given atom.
	 */
	private int next (int atomId) {
		
		int id = -1;
		
		if (this.atoms.get(atomId).getiEPSize() > 1){
		    this.atoms.get(atomId).dequeueFromiEP(); //removes and returns the first value in iep
		    id = this.atoms.get(atomId).getiEPValue(0);
		    
		} else {
			this.atoms.get(atomId).nullifyiEP();
			id = -1;
		}
		return id;
	}

	/*
	 * This function returns the next expected peak of the given atom
	 */
	private int getId(int atomId) {
	    int id;
        
	    id = this.atoms.get(atomId).getiEPValue(0);
        
		return id;
	}


	/*
	 * This function returns expected peak which could not be assigned on an earlier attempt
	 */
	private int getRejectedExpectedPeak() {
		double rand, threshold;
		int i, expId = -1, specId, expIdInSpec;
		
		rand = this.getRandomNumber();
		
		threshold = this.factorThreshold * rand;
		
		i = 0;
		reject: while (true){
			while (true){
				expId = -1;
				if (i == 100)
					break reject;
				expId = this.assignment.dequeueRejected();
				if (expId == -1) 
					break reject;
				if (this.assignment.getStatus(expId) == "PR"){
					break;
				}
			}	
			specId = this.expectedPeakSpecId[expId];
			expIdInSpec = this.expectedPeakExpIdInSpec[expId];
			
			
			if (this.allSpectra.get(specId).getExpectedPeak(expIdInSpec).getPeakVolumeFromFile() >= threshold ){
				break reject;
			}
			this.assignment.addRejectedPeak(expId);
			threshold = threshold * 0.95;
			i += 1;
		}   
		
		if ( i > 10) {
			this.factorThreshold = this.factorThreshold * 0.95;
		}else if (i < 3 && expId != -1) {
			this.factorThreshold =  this.factorThreshold * 1.05;
		}
		return expId;
	}

	/*
	 * This function makes all assignments that are possible for the unassigned expected peaks
	 */
	protected void makeAssignment(){
		int expectedPeakId; //id of the expected peak
		int measuredPeakId; //id of the measured peak
		
		while(unassignedExists(this.assignment.getUnassigned())){
			expectedPeakId = this.assignment.removeFromUnassigned();
			measuredPeakId = -1;
			
			if (this.assignmentType == 2) { 
				measuredPeakId = getPeakToOptimize(expectedPeakId);
			}else{ 
				measuredPeakId = getPeak(expectedPeakId);
			}
			
			if (measuredPeakId >= 0){
				assignExpPeakToMeasPeak(expectedPeakId,measuredPeakId);
				
			} else {  // measId = -1, so no measured peak is found.
				rejectPeakAssignment(expectedPeakId);
			}
		}
	
		this.assignment.resetUnassigned();

	}

	/*
	 * This function randomly mixes the existing values of the given array.
	 */
	private void randomizeArrayValues(int[] vals) {
		int begin, end, sizeP, ipeaks, changePos, temp;
		double random;		
		begin = 0;
		end = vals.length -1;
		sizeP = end - begin + 1;
		
		if (sizeP > 0) {
			for (ipeaks = begin; ipeaks <= end; ipeaks++) {
				random = this.getRandomNumber();
				changePos = (int) Math.floor(random* sizeP) + begin;
				temp = vals[ipeaks];
				vals[ipeaks] = vals[changePos];
				vals[changePos] = temp;		
			}			
		}
	}
	
	/*
	 * This function quicksorts the input and writes the sorted values into vals
	 */
	public void quickSort(int[] vals){
		int i;
		int[] position;
		int[] tempValue;
		
		int valsSize = vals.length;
		position = new int[valsSize];
		tempValue = new int[valsSize];
		
		for (i = 0; i < valsSize; i++) {
			position[i] = i;
			tempValue[i] = vals[i];
		}

		this.qs(vals, position, 0, valsSize-1);
		
		for (i = 0; i < valsSize; i++) {
			vals[i] = tempValue[position[i]];
		}
	}
	
	/*
	 * Quick sort algorithm. It sorts the values array and saves its indices in positions array.
	 */
	public void qs(int[] values, int[] positions, int beginIndex, int endIndex){
		
		int i, j, help;
		double medium;
		
	
		i = beginIndex;
		j = endIndex;
		medium = values[positions[(i+j)/2]];
		
		while (i <= j){
			while (values[positions[i]] < medium)
				i ++;
			while (medium < values[positions[j]])
				j--;
			if (i <= j){
				help = positions[i];
				positions[i] = positions[j];
				positions[j] = help;
				i++;
				j--;
			}
		}
		if (beginIndex < j)
			qs(values, positions, beginIndex, j);
		if (i < endIndex)
			qs(values, positions, i, endIndex);
	
	}
	
		/*
		 * This function is the quick sort algorithm. It sorts the values array and saves its indices in positions array.
		 */
		public void qs(double[] values, int[] positions, int beginIndex, int endIndex){
			
			int i, j, help;
			double med;
			i = beginIndex;
			j = endIndex;
			med = values[positions[(i+j)/2]];
			
			while (i <= j){
				while (values[positions[i]] < med)
					i ++;
				while (med < values[positions[j]])
					j--;
				if (i <= j){
					help = positions[i];
					positions[i] = positions[j];
					positions[j] = help;
					i++;
					j--;
				}
			}
			if (beginIndex < j)
				qs(values, positions, beginIndex, j);
			if (i < endIndex)
				qs(values, positions, i, endIndex);
		
			
		}
	
	
	/*
	 *	This function rejects the current assignment of the given expected peak. 
	 */
	private void rejectPeakAssignment(int expectedPeakID) {
		int spectrumId, expIdInSpec,dimension;
		spectrumId = this.expectedPeakSpecId[expectedPeakID];
		expIdInSpec = this.expectedPeakExpIdInSpec[expectedPeakID];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		
		this.assignment.addRejectedPeak(expectedPeakID);
		
		
	}
	/*
	*	This function checks the unassigned list and returns true if there is any value in this list. 
	*	It returns false otherwise. 
	*/
	public boolean unassignedExists(ArrayList<Integer> listUnassigned){
		
		if (listUnassigned.size() > 0)
			return true;
		else 
			return false;
	}

	/*
	 * This function assigns the expected peak to the measured peak that are given as input. 
	 */
	private void assignExpPeakToMeasPeak(int expPeakId, int measuredPeakId) {

		int spectrumId, expIdInSpec, iterator, dimension, atomsIndexInChemShift, atomIndexInSpec, i;
		boolean isAssigned;
		boolean [] adapt;
		int atomId, expectedPeakId;

		spectrumId = this.expectedPeakSpecId[expPeakId];
		expIdInSpec = this.expectedPeakExpIdInSpec[expPeakId];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		adapt = new boolean[dimension];
		

		for (i = 0; i < dimension; i++)
			adapt[i] = true;
		
		isAssigned = false;
		
		for (iterator = 0; iterator < dimension; iterator++){
			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(iterator).getIndex();
			this.atoms.get(atomsIndexInChemShift).setAssigned(true);
			this.allSpectra.get(spectrumId).getAtoms().get(atomIndexInSpec).setAssigned(true);
			atomId = atomsIndexInChemShift;
			if (adapt[iterator] && this.allSpectra.get(spectrumId).getAtom(atomIndexInSpec).getMeasuredPeakAvgSum() == 0){
				
				startPeakId(atomId);
				
				for (i = 0; i < this.atoms.get(atomsIndexInChemShift).getExpectedPeakSize(); i++){ //for all expected peaks in landscape for this atom
					expectedPeakId = this.atoms.get(atomsIndexInChemShift).getExpectedPeak(i).getAssignmentIndex();
					if ( expPeakId != expectedPeakId && ((this.assignment.getAssignment(expectedPeakId) == -1)&&(this.assignment.getRejectNumber(expectedPeakId) > 0))){
						this.assignment.addToUnassigned(expectedPeakId);
						this.assignment.setExpStatusToUnassigned(expectedPeakId);			
					}
				}
			}
		} 
		
		setPeakAssignment(expPeakId, measuredPeakId);
		
	}


	/*
	 * This function assigns the given measured peak to the given expected peak.
	 */
	private void setPeakAssignment(int expectedPeakId, int measuredPeakID) {
		
		int specId, expIdInSpec, dimension, idim;
		int atomIndexInSpec, atomsIndexInChemShift;
		specId = this.expectedPeakSpecId[expectedPeakId];
		expIdInSpec = this.expectedPeakExpIdInSpec[expectedPeakId];
		dimension = this.allSpectra.get(specId).getDimension();
		
		
		for (idim = 0; idim < dimension; idim ++){
			atomsIndexInChemShift =  this.allSpectra.get(specId).getExpectedPeaks().get(expIdInSpec).getAtom(idim).getAtomsIndex();
			atomIndexInSpec = this.allSpectra.get(specId).getExpectedPeaks().get(expIdInSpec).getAtom(idim).getIndex();
			
           
			addFrequency(expectedPeakId, measuredPeakID, this.allSpectra.get(specId).getChemicalShiftValues(idim,measuredPeakID), idim);
		}
	
		setMeasuredPeakIdInExpectedPeak(expectedPeakId, measuredPeakID);

	}


	
	/*
	 * This function writes measured peak ID  to expected peak
	 */
	private void setMeasuredPeakIdInExpectedPeak(int expectedPeakId, int measuredPeakId) {
		
		int spectrumId, expectedPeakIdInSpectrum, dimension;

		spectrumId = this.expectedPeakSpecId[expectedPeakId];
		expectedPeakIdInSpectrum = this.expectedPeakExpIdInSpec[expectedPeakId];
		dimension = this.allSpectra.get(spectrumId).getDimension();
		this.assignment.setAssignment(expectedPeakId, measuredPeakId);

		this.allSpectra.get(spectrumId).increaseAssignedPeakCount(measuredPeakId);
		
		
	}

	/*
	 * This function updates the chemical shift frequencies of the given expected peak with the values 
	 * from the given measured peak.
	 * This function is called after an assignment done, to update all atoms' chemical shift values with this assignment.
	 */
	private void addFrequency(int expectedPeakId, int measuredPeakId, double chemShift, int dimension) {
		int spectrumId, expIdInSpec, atomsIndexInChemShift, atomIndexInSpec;
		spectrumId = this.expectedPeakSpecId[expectedPeakId];
		expIdInSpec = this.expectedPeakExpIdInSpec[expectedPeakId];
		atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(dimension).getAtomsIndex();
		atomIndexInSpec = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(dimension).getIndex();
		this.addFrequencyCounter = this.addFrequencyCounter +1;
		
		//add values to the landscape atom and set atom to assigned
		this.atoms.get(atomsIndexInChemShift).addFrequencyToChemShiftAtom(chemShift);
		this.atoms.get(atomsIndexInChemShift).setAssigned(true);
		
		//add values to the spectrum atom
		this.allSpectra.get(spectrumId).getAtom(atomIndexInSpec).addFrequencyToSpectrumAtom(chemShift);
		this.allSpectra.get(spectrumId).getAtom(atomIndexInSpec).setAssigned(true);
	}


	/*
	 * This function finds and returns the measured peak, that fits best to expected peak,expIdGlobal
	 * If no measured peak is found, then -1 is sent. 
	 */
	private int getPeak(int expIdGlobal){
		int peakId = -1;
		int expId; //expected peak Id in dedicated spectrum.expectedPeaks
		int spectrumId; //spectrumId that the expected peak belongs to.
		int expSize; // number of expected peaks of the specific spectrum
		int dimension; //dimension of the specific spectrum (#atoms in each peak)
		int i, j,atomsIndex;
		Atom at = new Atom();
		double expectedFrequency, expectedFrequencyDeviation, chemShift, deviation;
		spectrumId = this.expectedPeakSpecId[expIdGlobal];
		expId =  this.expectedPeakExpIdInSpec[expIdGlobal];
		
		expSize = this.allSpectra.get(spectrumId).getExpectedPeaksSize(); //exp peak size of that spectrum
		dimension = this.allSpectra.get(spectrumId).getDimension();

		for (i = 0; i < dimension; i++){
			getFrequencyRange (spectrumId, expId, i, false, false);
			getExpectedFrequency(spectrumId, expId, i);
			getMeasuredFrequency(spectrumId, expId, i);

		}
		boolean [] hold = new boolean[dimension];
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expId).getAtom(i).getAtomsIndex();
			hold[i] = this.atoms.get(atomsIndex).getAssigned(this.assignment.getAssignmentWhole());
		}
		peakId = selectPeak(spectrumId, expId, hold);
	
		return peakId;
		
	}
	
	/*
	 * This function selects and returns the ID of the expected peak whose assignment to be improved. 
	 */
	private int getPeakToOptimize(int expIdGlobal){
		int peakId = -1;
		int expectedPeakId; //expected peak Id in dedicated spectrum.expectedPeaks
		int spectrumId; //spectrumId that the expected peak belongs to.
		int expectedPeakCounter; // number of expected peaks of the specific spectrum
		int dimension; //dimension of the specific spectrum (#atoms in each peak)
		int i, j,k,atomsIndex,atomsIndexInChemShift;
		Atom at = new Atom();
		double random;
		boolean [] assignment;

		spectrumId = this.expectedPeakSpecId[expIdGlobal];
		expectedPeakId =  this.expectedPeakExpIdInSpec[expIdGlobal];
		this.allSpectra.get(spectrumId).getExpectedPeak(expectedPeakId).setAssignmentIndex(expIdGlobal);
		dimension = this.allSpectra.get(spectrumId).getDimension();
		
		
		expectedPeakCounter = this.allSpectra.get(spectrumId).getExpectedPeaksSize(); //exp peak size of that spectrum
		assignment = new boolean[dimension];
		
		for (i = 0; i < dimension; i++){

			atomsIndexInChemShift =  this.allSpectra.get(spectrumId).getExpectedPeaks().get(expectedPeakId).getAtom(i).getAtomsIndex();
			assignment[i] = this.atoms.get(atomsIndexInChemShift).getAssigned(this.assignment.getAssignmentWhole());
			getFrequencyRange (spectrumId, expectedPeakId, i, false, false);
			getExpectedFrequency(spectrumId, expectedPeakId, i);
			getMeasuredFrequency(spectrumId, expectedPeakId, i);
	
		}
		
		if (countTrues(assignment) == dimension) {
			peakId = selectPeak(spectrumId, expectedPeakId, assignment);
   		 	
		} else {
			peakId = selectedEvolutionary(spectrumId, expectedPeakId, assignment, assignment);
			if (peakId == -1 && this.temp > 0) {
				random = getRandomNumber();
				random = Math.pow(random, this.temp);
				if (random <= 0.5) {
					if (this.allSpectra.get(spectrumId).getExpectedPeak(expectedPeakId).equalsAllocated()) {
						peakId = getEqualPeakId(spectrumId, expectedPeakId, assignment);
					}

					boolean assigned = false;
					boolean diagonal = false;
					for (k = 0; k < dimension; k++) {
						assigned = assigned || assignment[k];
						diagonal = diagonal ||  this.allSpectra.get(spectrumId).getExpectedPeak(expectedPeakId).getDiagonal()[k];
					}
					
	                if (peakId == -1 && (! (!assigned && diagonal)) && random <= 0.3) {
	                	peakId = selectPeak(spectrumId, expectedPeakId, assignment);
	                }
				}
			}
		}
		return peakId;
		
	}
	/*
	 * This function returns the assigned measured peak id for the given expected peak id in the given individual. 
	 */
	private int getMeasIdForExpIdInIndividual(int individualId, int expId) {
		
		int measuredPeakId = -1;
		if (expId != -1) {
			measuredPeakId = Integer.parseInt(this.population.get(individualId).getVariable(expId).toString());
		}
		
 		return measuredPeakId;
	}
	
	/*
	 * This function selects the measured peak for the given expected peak as per their probability.
	 */
	private int selectedEvolutionary(int spectrumId, int expectedPeakId, boolean[] assignment, boolean[] hold) {
	
		int peakId, individual, value = 0, degeneracy = 0, activeLength, assignmentId, iterator, assignedPeakCounter;
		int i, dimension, atomsIndex, expIdInSpec, expIdGlobal, assignmentIterator, individualId;
		double probability, chemShift, random;
		ArrayList<Integer> peakValuesList = new ArrayList<Integer>();
		ArrayList<Integer> degeneracyPeakList = new ArrayList<Integer>();
		ArrayList<Double> probabilityDegeneracy = new ArrayList<Double>();
		ArrayList<Double> selectionProbability = new ArrayList<Double>();
		double [][] searchSpace ;
		double[] down, high;
		boolean goOn = false;
		activeLength = this.populationActive.length;
		
		expIdInSpec = expectedPeakId;
		expIdGlobal = this.allSpectra.get(spectrumId).getExpectedPeak(expectedPeakId).getAssignmentIndex();
		dimension = this.allSpectra.get(spectrumId).getDimension();
		searchSpace = new double[dimension][2];
		down = new double[dimension];
		high = new double[dimension];
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
			searchSpace[i]= this.atoms.get(atomsIndex).getSearchSpace();
		}
		peakId = -1;
		individual = -1;
		
		peak: for (assignmentIterator = 0; assignmentIterator < activeLength; assignmentIterator++) {
			individualId = this.populationActive[assignmentIterator];
		
			assignmentId = getMeasIdForExpIdInIndividual(individualId, expIdGlobal);
			if (assignmentId != -1) {
				for (i = 0; i < dimension; i++) {
					down[i] = searchSpace[i][0];
					high[i] = searchSpace[i][1];
				}
				
				for (iterator = 0; iterator < dimension; iterator++) {
					chemShift = this.allSpectra.get(spectrumId).getChemicalShiftValues(iterator, assignmentId);
					goOn = true;
					if (allIsTrue(assignment) && (chemShift < down[iterator] || chemShift > high[iterator]))
						goOn = false;
					
				}
				if (goOn) {
					
					assignedPeakCounter = this.allSpectra.get(spectrumId).getAssignedPeakCount(assignmentId);
					if (assignedPeakCounter == 0) {
						probability = calculateProbability(spectrumId, assignmentId, expIdInSpec, populationActive[assignmentIterator]);
						peakValuesList.add(value,assignmentId);
       					selectionProbability.add(value,probability);
       					value += 1;
       					
				
					} else if (value == 0) {
						probability = calculateProbability(spectrumId, assignmentId, expIdInSpec, populationActive[assignmentIterator]);
						degeneracyPeakList.add(degeneracy, assignmentId);
						probabilityDegeneracy.add(degeneracy,probability);
						degeneracy += 1;
					}
				}
			}
		}
		
		 if (value > 0){
			individual =  select (selectionProbability);
				
		    peakId = peakValuesList.get(individual);
		    
		 } else if (degeneracy > 0  ) {  
			 random = this.getRandomNumber();
			 random = random * this.probabilityDegeneracy;
			 individual = select (probabilityDegeneracy);
				if ((random < this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).getPeakVolumeFromFile()&& !(false))||
						(allIsTrue(hold))) {
					individual = select (probabilityDegeneracy);
					peakId = degeneracyPeakList.get(individual);
					
				}
			}
		return peakId;
	}

	/*
	 * This function returns true if values of the given array is true.
	 * It returns false otherwise.
	 */
	private boolean allIsTrue(boolean[] array) {
		for (int i = 0; i < array.length; i ++) {
			if (array[i]== false)
				return false;
		}
		return true;
	}

	/*
	 * This function seeks an expected peak that is equal to the input expected peak. 
	 */
	private int getEqualPeakId(int spectrumId, int expectedPeakId, boolean[]assignment) {
		
		int activeLength = this.populationActive.length;
		int i , expIdInSpec, peakId = -1, atomsIndex, pExpPeak;
		double [][] searchSpace;
		double [] down, high;
		int dimension = this.allSpectra.get(spectrumId).getDimension();
	
		searchSpace = new double[dimension][2];
		down = new double[dimension];
		high = new double[dimension];
		expIdInSpec = expectedPeakId;
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
			searchSpace[i]= this.atoms.get(atomsIndex).getSearchSpace();
		}
		pExpPeak = expectedPeakId;
			
		if (this.allSpectra.get(spectrumId).getExpectedPeak(pExpPeak).equalsAllocated()){
			for (i = 0; i < dimension; i++) {
				down[i] = searchSpace[i][0];
				high[i] = searchSpace[i][1];
			}
			if (this.allSpectra.get(spectrumId).getExpectedPeak(pExpPeak).getEquivPeakSize() == 0 ) {
	
			}
		
		}
		
		return peakId;
	}
	
	/*
	 * This function returns the current search space. 
	 */
	private double getRange(int spectrumId, int peakId) {
	
		double range = 0;
		int residueFirst, residueSecond, atomIdFirst, atomIdSecond;
		
		for (atomIdFirst = 0; atomIdFirst < this.allSpectra.get(spectrumId).getExpectedPeak(peakId).getAtoms().length; atomIdFirst++) {
			for (atomIdSecond = 0; atomIdSecond < this.allSpectra.get(spectrumId).getExpectedPeak(peakId).getAtoms().length; atomIdSecond++) {
				residueFirst = this.allSpectra.get(spectrumId).getExpectedPeak(peakId).getAtom(atomIdFirst).getResidueNumber();
				residueSecond = this.allSpectra.get(spectrumId).getExpectedPeak(peakId).getAtom(atomIdSecond).getResidueNumber();
				range = Math.max(range, Math.abs(residueFirst - residueSecond));
			}
		}
		return range;
	}
	
	/*
	 * This function searches the search space and selects the measured peak to be assgined for the given expected peak.
	 */
	private int selectPeak(int spectrumId, int expectedPeakId, boolean[] save) {
		int peakId, in, value = 0, degree = 0, peakNumber, iterator;
		int i, dimension, atomsIndex, index, expIdInSpec;
		double probability;
		ArrayList<Integer> peakIndexList = new ArrayList<Integer>();
		ArrayList<Integer> degeneracyList = new ArrayList<Integer>();
		ArrayList<Double> probabilityList = new ArrayList<Double>();
		ArrayList<Double> probabilitySelectionList = new ArrayList<Double>();
		double [][] searchSpace ;
		ArrayList<Integer> peakInList = new ArrayList<Integer>();

		peakNumber = 0;		
		expIdInSpec = expectedPeakId;
		
		dimension = this.allSpectra.get(spectrumId).getDimension();
		
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
		}
		
		searchSpace = new double[dimension][2];
		
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
			searchSpace[i]= this.atoms.get(atomsIndex).getSearchSpace();
		}
		peakId = -1;
		in = -1;

		if (countTrues(save) > dimension - 2){			
			peakInList.clear();
			peakInList = traverseSearchSpace(spectrumId, expIdInSpec, searchSpace);
			if (peakInList.size() > 0){
				value = 0;
				degree = 0;
				peakNumber = peakInList.size();
				peakIndexList = new ArrayList<Integer>(peakNumber);
				degeneracyList = new ArrayList<Integer>(peakNumber);
				probabilitySelectionList = new ArrayList<Double>(peakNumber);
				probabilityList = new ArrayList<Double>(peakNumber);
			}
			
			for (iterator = 0; iterator < peakNumber; iterator++){
				
				if (this.allSpectra.get(spectrumId).getAssignedPeakCount(peakInList.get(iterator)) == 0){
					probability = calculateSelectionProbability(spectrumId, peakInList.get(iterator), expIdInSpec);
					peakIndexList.add(value,peakInList.get(iterator));
					probabilitySelectionList.add(value,probability);
					value += 1;
	                
				} else if (value == 0){
					 probability = calculateSelectionProbability(spectrumId, peakInList.get(iterator), expIdInSpec);
					 degeneracyList.add(degree, peakInList.get(iterator));
					 probabilityList.add(degree,probability);
					 degree += 1;
				}
			}
			
			if (value > 0){
		        in =  select (probabilitySelectionList);
                peakId = peakIndexList.get(in);

			} else if (degree > 0){
				in = select (probabilityList);
				peakId = degeneracyList.get(in);
			}
			
		} else  {
			peakId = search(spectrumId, expIdInSpec, save);
			
		}
		
		return peakId;
	}

	/*
	 * This function selects the measured peak depending on the given probability list. 
	 */

	private int select(ArrayList<Double> probabilitiesInputList) {
		
		int index, inputSize, start, last;
		double[] probabilitySum;
		double probValue, random, difference;
		
		inputSize = probabilitiesInputList.size();
		
		if (inputSize == 1)
			index = 0;
		else{
			probabilitySum = new double[inputSize];
			probabilitySum[0] = probabilitiesInputList.get(0);
			for (int i = 1; i < inputSize; i++){
				probabilitySum[i] = probabilitySum[i-1] + probabilitiesInputList.get(i);
			}
			
			probValue = probabilitySum[inputSize-1];
			random = this.getRandomNumber();
			
			random = random * probValue;
			index = 0;
			start = 0;
			last = inputSize-1;
			difference = 1;
			while(start <=last && difference != 0){
				index = (start + last)/ 2;
				difference = probabilitySum[index] - random;
				if (difference < 0)
					start = index +1;
				else if (difference < probabilitiesInputList.get(index))
					difference = 0;
				else
					last = index -1;
			}
		}
	    return index;
	}
	
	/*
	 * This function calculates the probability for the given expected peak.
	 */
	private double calculateProbability(int spectrumId, int assignmentId, int expIdInSpec, int position){
		
	   double scale, scaleValue, probability;
	   int residueId;
	   
	   residueId = this.allSpectra.get(spectrumId).getExpectedPeak(expIdInSpec).getResidueNumber();
	   residueId = residueId - 1; //first residue Id is 1, so 1 is subtracted
	   
	   scale = Math.abs(this.maximum[residueId]- this.minimum[residueId]);
	    
	   scaleValue = this.middleScore[position][residueId] - this.minimum[residueId];
	    
	   if (scale == 0 )
		   probability = 1;
	   else
		   probability = scaleValue / scale;

	   probability = Math.pow(probability, 10);
	   probability = Math.max(0.0001, probability);
	   
	   return probability;
	}
	
	/*
	 * This function calculates the selection probability of the given measured peak and the expected peak.
	 */
	private double calculateSelectionProbability(int spectrumId, int measuredPeakId, int expectedPeakId) {
		
		int dimension,  atomsIndex,  iterator, position;
		double probabilityTotal = 1, chemShift, deviationMeasuredFreq, measuredFreq, probability, expectedChemShiftDev, expectedChemShift;
		
		dimension = this.allSpectra.get(spectrumId).getDimension();
		
		for (iterator = 0; iterator < dimension; iterator++){
			chemShift = this.allSpectra.get(spectrumId).getChemicalShiftValues(iterator, measuredPeakId);
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expectedPeakId).getAtom(iterator).getAtomsIndex();
			measuredFreq = this.atoms.get(atomsIndex).getMeasFreqDynamic();
			if (measuredFreq != this.NNN){
				
				deviationMeasuredFreq = this.atoms.get(atomsIndex).getMeasDevDynamic();
				position = (int) Math.floor(Math.abs((measuredFreq- chemShift)/deviationMeasuredFreq)*this.residueAssignment);
				if (position > this.residueAssignment){
					position = this.residueAssignment;
				}
				probability = this.assignmentDistance[position]; 
				
			} else{
				
				expectedChemShiftDev = this.atoms.get(atomsIndex).getDynamicExpectedFreqDeviation();
				expectedChemShift = this.atoms.get(atomsIndex).getDynamicExpectedFrequency();
				if (expectedChemShiftDev != 0){
					position = (int) Math.floor(Math.abs((expectedChemShift - chemShift)/expectedChemShiftDev /4)* this.residueExpression);
					if (position > this.residueExpression){
						position = this.residueExpression;
					}
				} else {
					position = this.residueExpression;
				}
				probability = this.expScore[position];
				      
			}
		    
			probabilityTotal = probabilityTotal * probability; 
	
		}
		
		probabilityTotal = probabilityTotal * 1; 
		
		return probabilityTotal;
	}

	/*
	 *  This function searches the current search space to find a suitable assignmetn for the given expected peak 
	 *  during the population initialization.
	 */
	private int search (int spectrumId, int expectedPeakId, boolean[] save){
		int peak = -1, i, atomsIndex, peakId, ipeak;
		int is = 0, dimension, iterator, total;
		int [] freePeaks;
		double [][] searchSpace;
		double random, deviation;

		double [][] ranges;
		ArrayList<Integer> positionList = new ArrayList<Integer>();
		
		positionList.clear();
		
		dimension = this.allSpectra.get(spectrumId).getDimension();
		ranges = new double[dimension][2];
		searchSpace = new double [dimension][2];
		
				
		for (i = 0; i < dimension; i++){
			atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expectedPeakId).getAtom(i).getAtomsIndex();
			ranges[i]= this.atoms.get(atomsIndex).getSearchSpace();
		}
		
		peakId = 0;
		for (is = -1; is <= 1; is++){
			for (iterator = 0; iterator < dimension; iterator++){
				atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expectedPeakId).getAtom(iterator).getAtomsIndex();
				if (save[iterator] == true){
					searchSpace[iterator] = this.atoms.get(atomsIndex).getSearchSpace();					
				} else {
					random = getRandomValue(atomsIndex);
					deviation = this.atoms.get(atomsIndex).getDynamicExpectedFreqDeviation();
					
					searchSpace[iterator][0] = random - deviation * Math.pow(2.0, is); 
					searchSpace[iterator][1] = random + deviation * Math.pow(2.0, is); 
					
					if (searchSpace[iterator][0] < ranges[iterator][0] || searchSpace[iterator][1]  > ranges[iterator][1]){
						searchSpace[iterator] = ranges[iterator];
					}
				}
			}
			
			positionList.clear();
			positionList = traverseSearchSpace(spectrumId, expectedPeakId, searchSpace);
			
			peakId = sumUnassignedPeaks(spectrumId, positionList);
			if (peakId > 0){
				break;
			}
		}
		
		if (peakId == 0){
			searchSpace = ranges; 
			positionList.clear();
			positionList = traverseSearchSpace(spectrumId, expectedPeakId, searchSpace);
			peakId = sumUnassignedPeaks(spectrumId, positionList);
		}
		
		if (positionList.size()> 0 && peakId > 0){
			total = 0;
			freePeaks = new int[peakId];
			for (ipeak= 0 ; ipeak < positionList.size(); ipeak++){
				if (this.allSpectra.get(spectrumId).getAssignedPeakCount(positionList.get(ipeak)) == 0){
					freePeaks[total] = positionList.get(ipeak);
					total ++;
				}
			}
			random = this.getRandomNumber();
			peak = freePeaks[(int) Math.floor(peakId*random)];
		}

		return peak;
	}
	
	/*
	 * This function returns the total number of unassigned peaks in the given spectrum.
	 */
	private int sumUnassignedPeaks( int spectrumId, ArrayList<Integer> positionList){
		int unassignedPeakCounter = 0;
		int i;
		
		
		if (positionList.size() > 0){
			for (i = 0; i < positionList.size(); i++){
				if (this.allSpectra.get(spectrumId).getAssignedPeakCount(positionList.get(i))== 0)
					unassignedPeakCounter ++;
			}
		}
		
		return unassignedPeakCounter;
	}
	
	/*
	 * This function traverses the search space to func the possible positions for an assignment. 
	 */
	private ArrayList<Integer> traverseSearchSpace(int spectrumId, int expectedId, double[][] ranges) {
		int i;
		ArrayList<Integer> positions = new ArrayList<Integer>();
		int[] positionRight = new int[this.allSpectra.get(spectrumId).getMeasuredPeaksSize() + 1];
		positions.clear();
		
		positionRight[0] = 1;
		positionRight = getRange(this.allSpectra.get(spectrumId).getTree(), ranges, 0, positionRight);
		if (positionRight[0] > 1){ //if there is already value assigned  (except the first one, since it shows only the next empty index, so 
			for (i = 1; i< positionRight[0]; i ++){
				positions.add(positionRight[i]);
			}
		}
		
		return positions;
	
	}
	
	/*
	 * This function returns the current search range. 
	 */
	private int[] getRange(Tree myTree, double[][] ranges, int dimension, int[] positions){
		int rangeLength, sNode, left, right, root, start, end, childRight, childLeft;
		double leftBord, rightBord, valueNode;
		
		int [] result = new int[3];
		
		
		
		rangeLength = ranges.length;
		leftBord = ranges[dimension][0];
		rightBord = ranges[dimension][1];
		

		result = splitFind(myTree, leftBord, rightBord);
		sNode = result[0];
		left = result [1];
		right = result [2];
		
		
		if (left == right){
			positions = checkChild(myTree, sNode, dimension, rangeLength, positions, ranges);
		} else if ((myTree.getMeasuredPeakInSpectrum(left) >= leftBord) &&  (myTree.getMeasuredPeakInSpectrum(right)<= rightBord)){
			if ((dimension +1) < rangeLength){
				getRange(myTree.getFollowingDimension(sNode), ranges, (dimension+1), positions);
			} else{
				positions = findPosition(myTree, positions, left, right);
			}
		} else {
			start = left;
			end = sNode; 
			while (start < end){
				root = (start + end) / 2;
				valueNode = myTree.getMeasuredPeakInSpectrum(root);
				if (leftBord <= valueNode){
					if ((dimension+1) == rangeLength){
						positions = findPosition(myTree, positions, (root +1), end);
					} else{
						if ((root +1) == end){
							positions = checkChild(myTree.getFurtherDimension(end), 0, dimension+1, rangeLength, positions, ranges);
						} else{
							childRight = (root +1+end)/2;
							positions = getRange(myTree.getFollowingDimension(childRight), ranges, dimension+1, positions);
						
						}
					}
					end = root;
					
				} else{
					start = root +1;
				}
			}
			root = start;
			positions = checkChild(myTree, root, dimension, rangeLength, positions, ranges);
			start = sNode + 1;
			end = right;
			root = start;
			while (start < end){
				root = (start + end)/2;
				valueNode = myTree.getMeasuredPeakInSpectrum(root);
				if (rightBord >= valueNode){
					if ((dimension+1) == rangeLength){
						positions = findPosition(myTree, positions, start, root);
					} else {
						if (root == start){
							positions = checkChild(myTree.getFurtherDimension(root), 0, dimension+1, rangeLength, positions, ranges);
						} else{
							childLeft = (root + start)/ 2;
							getRange(myTree.getFollowingDimension(childLeft), ranges, dimension+1, positions);
						}
					}
					start = root +1;
				} else {
					end = root;
				}
			}
			root = start;
			positions = checkChild(myTree, root, dimension, rangeLength, positions, ranges);
		}
			
		
		return positions;
	}

	/*
	 * This function returns the index of the given leaf value.
	 */

	private int[] checkChild(Tree myTree, int leaf, int dimension, int rangeLength, int[] positions, double[][] ranges) {
		
		int d, l;
		
		d = dimension;
		Tree pTree = new Tree(myTree);
		l = leaf;
	
		while (true){
			if (!((pTree.getMeasuredPeakInSpectrum(l) <= ranges[d][1])&&(pTree.getMeasuredPeakInSpectrum(l) >= ranges[d][0]))){
				return positions;
			}
			if ((d+1) == rangeLength){
				break;
			}
			 d += 1;
			
			 pTree = pTree.getFurtherDimension(l);
				
			 l = 0;
		}
		positions = findPosition( pTree, positions, l, l);

		return positions;
	}

	/*
	 * This function returns the positions of the given indices in the tree.
	 */
	private int[] findPosition(Tree myTree, int[] positions, int first, int last) {
		int numPoints;
		int i, initial;
		
		numPoints = last - first +1;
		initial = positions[0];
		
		for (i = 0; i < numPoints; i++){
			positions[i+initial] = myTree.getOrientationOf(first+i); 
		}
		positions[0] = positions[0] + numPoints;

		
		return positions;
	}


	/*
	 * This function searches for the minimum and maximum values in the tree
	 */
	private int[] splitFind(Tree myTree, double minimumValue, double maximumValue){
		int []result = new int[3];
		int node, first, last, root;
		double value;
		
		node = -1;
		
		first = 0;
		last = myTree.getMeasuredPeaks().size() -1; //if 14 elements, then array[13] is the last value since it starts from 0
		
		while (first < last){
			root = (first + last)/2;
			value = myTree.getMeasuredPeakInSpectrum(root);
			if (minimumValue > value){
				first = root +1;
			} else if (maximumValue < value){
				last = root;
			}else{
				node = root;
				result [0] = node;
				result [1] = first;
				result [2] = last;
				return result;
			}
		}
		
		node = first;
		result [0] = node;
		result [1] = first;
		result [2] = last;
		
		return result;
		
	}
	
	private double getRandomValue(int atomsIndex){
		double result;
		
		result =  getRandomDeviation(atomsIndex);
		
		return result;
	}
	
	private double getRandomDeviation(int atomsIndex){
		double get_randValDev, randVal;
		double val, dev;
		
		val = this.atoms.get(atomsIndex).getDynamicExpectedFrequency();
		dev = this.atoms.get(atomsIndex).getDynamicExpectedFreqDeviation();
		randVal = normalDistributionRandomValue();
		get_randValDev = val + dev * randVal;
		return get_randValDev;
		
	}
	
	/*
	 * This function searches for the assigned chemical shift value and its accuracy for the given expected peak.
	 */
	private void getMeasuredFrequency(int spectrumId, int expectedId, int i) {
		int atomsIndex, index;
		double frequency = NNN;
		double accuracy = 0;
		double[] result = new double[2];
		int expIdInSpec;

		expIdInSpec = expectedId; //this.expectedPeakExpIdInSpec[expId];
		atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
		index = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getIndex();
		this.atoms.get(atomsIndex).getCsTableFrequency();
		this.allSpectra.get(spectrumId).getAtoms().get(index);
		
		
		if (this.allSpectra.get(spectrumId).getAtoms().get(index).getAssigned(this.assignment.getAssignmentWhole()) == true){ 
			frequency = this.allSpectra.get(spectrumId).getAtoms().get(index).getMeasuredPeakAvgSum();
			accuracy = this.allSpectra.get(spectrumId).getAtoms().get(index).getAccuracy();
		
		} else if (this.atoms.get(atomsIndex).getAssigned(this.assignment.getAssignmentWhole())==true){
			frequency = this.atoms.get(atomsIndex).getMeasuredPeakAverage();
			accuracy = this.allSpectra.get(spectrumId).getAtoms().get(index).getAccuracy();
			
		}else{
			frequency = this.NNN;
			accuracy = 0;
		}
		
		result[0] = frequency;
		result[1] = accuracy;
		
		this.atoms.get(atomsIndex).setDynamicMeasuredFrequency(accuracy, frequency);
		this.allSpectra.get(spectrumId).getAtoms().get(index).setDynamicMeasuredFrequency(accuracy, frequency);
			
	}

	/*
	 * This function returns the total number of true values in the given array.
	 */
	private int countTrues(boolean[] array){
		
		int counter = 0;
		int i;
		
		for (i = 0; i < array.length; i++){
			if (array[i] == true)
				counter++;
		}
		return counter;
		
	}
	
	/*
	 * This function returns the expected chemical shift frequency from the CS table of the given expected peak.
	 */
	private void getExpectedFrequency(int spectrumId, int expectedPeakId, int i) {
		int atomsIndex, index;
		double expFreqVal, expFreqDev;
		double[] result = new double[2];
		int expIdInSpec;
		
		expIdInSpec = expectedPeakId; 
		atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
		index = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getIndex();
		
		expFreqDev = this.atoms.get(atomsIndex).getCsTableFrequency();
		expFreqVal = this.atoms.get(atomsIndex).getReferenceFrequency();
		result[0] = expFreqDev;
		result[1] = expFreqVal;
		
		this.atoms.get(atomsIndex).setDynamicExpectedFrequency(expFreqDev, expFreqVal);
		this.allSpectra.get(spectrumId).getAtoms().get(index).setDynamicExpectedFrequency(expFreqDev, expFreqVal);
		
					
	}


	/*
	 * This function assigns the lower and upper boundary of the expected chemical shift frequency of the i. atom fo the given expected peak.
	 */
	private void getFrequencyRange(int spectrumId, int  expectedPeakId, int i, boolean newStateValid, boolean newstate) {
		
		int atomsIndex; //index of the atom in chemicalShiftAssignment.atoms
		int index; //index of the atom in Spectrum.atoms
		double accuracy, randomDeviation, expectedChemShift, csTableFrequency, low = 0, up = 0, delete;
		boolean fixFreq = false;
		boolean state;
		int expIdInSpec;
		

		expIdInSpec = expectedPeakId;
		atomsIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getAtomsIndex();
		index = this.allSpectra.get(spectrumId).getExpectedPeaks().get(expIdInSpec).getAtom(i).getIndex();
		
	
		
		state = this.atoms.get(atomsIndex).getAssigned(this.assignment.getAssignmentWhole());
		if (newStateValid){
			state = newstate;
		}
		if (state){ //if atom is assigned in the landscape 
			
			if (this.allSpectra.get(spectrumId).getAtoms().get(index).getAssigned(this.assignment.getAssignmentWhole())){  //if atom is assigned in the same spectrum
				accuracy = this.allSpectra.get(spectrumId).getAtoms().get(index).getAccuracy();
				
                low = this.allSpectra.get(spectrumId).getAtoms().get(index).getMeasuredPeakAvgSum()- accuracy;
                up =  this.allSpectra.get(spectrumId).getAtoms().get(index).getMeasuredPeakAvgSum() + accuracy;
			} else{  //atom is not assigned in the same spectrum
				accuracy = this.allSpectra.get(spectrumId).getAtoms().get(index).getAccuracy();
				low = this.atoms.get(atomsIndex).getMeasuredPeakAverage() - accuracy;
		        up =  this.atoms.get(atomsIndex).getMeasuredPeakAverage() + accuracy;
			}
			
		} else{  //atom is not assigned in the landscape
			
			if (fixFreq){

			} else{ //if fix atom assignment is not given
				
				csTableFrequency = this.atoms.get(atomsIndex).getCsTableFrequency();
				expectedChemShift = this.atoms.get(atomsIndex).getReferenceFrequency();
				randomDeviation = getRandomDeviation(expectedChemShift,csTableFrequency);
				delete = factor * csTableFrequency + randomDeviation;
				low = expectedChemShift - delete;
				up = expectedChemShift + delete;
				
				
			}
			
		}
		
		this.atoms.get(atomsIndex).setSearchBoundary(low, up);
		this.allSpectra.get(spectrumId).getAtoms().get(index).setSearchBoundary(low, up);
	}
	
	/*
	 * This function returns the random deviation 
	 */
	private double getRandomDeviation (double expectedFrequency, double expectedFrequencyDeviation){
		double result, value;

		value = normalDistributionRandomValue();
		
		result = Math.abs(expectedFrequencyDeviation * value);
		
		return result;
	}
	

	
	/*
	 * This function returns random value according to normal distribution
	 */
	private double normalDistributionRandomValue(){
		double gaus;
		double first = 0, second = 0;
		double result = 0 ;

		if (this.gaussianStored){
			gaus = this.gaussian; 
			this.gaussianStored  = false;
		
		} else {
			while (true){
				first = this.getRandomNumber();
		
				second = this.getRandomNumber();
		
				first = 2*first - 1;
				second = 2*second - 1;
							
				
				result = Math.pow(first, 2) + Math.pow(second, 2);
				if (result > 0 && result < 1)
					break;
			}
			result = Math.sqrt(-2 * Math.log(result)/result);
			gaus = first * result;
			this.gaussian = second * result;
			this.gaussianStored = true;
		}
		return gaus;
		
	}
	

	/*
	 * This function reset all assignment related attributes in the landscape
	 */
	private void resetAssignment(){
		int i, j, expPeakCounter;
		resetAllPeakAssignments(); 
		
		this.assignment.reset();  // the new class is reset and unassigned list is randomized here
		resetAtomAssignments();
		
		expPeakCounter = 0;
		for (i = 0; i < this.allSpectra.size(); i++){  //for each spectrum
			for (j = 0; j < this.allSpectra.get(i).getAtoms().size(); j++){ //for each atom of that spectrum
				this.allSpectra.get(i).getAtoms().get(j).clearMeasuredFrequency();
				this.allSpectra.get(i).getAtoms().get(j).calculateAverageFrequency();
			}
			
			this.allSpectra.get(i).setInitialValues();
			this.allSpectra.get(i).setAssignmentIndexFirstExpPeak(expPeakCounter);
			this.allSpectra.get(i).setAssignmentIndexLastExpPeak(expPeakCounter+this.allSpectra.get(i).getExpectedPeaksSize()-1 );
			expPeakCounter += this.allSpectra.get(i).getExpectedPeaksSize();
		}
		
		
	}
	

	/*
	 * This function resets the measured peak frequency values in chemicalShift.atoms
	 */
	private void resetAllPeakAssignments() {
		int pExpPeak;
		int i, j, k, specId;
		
		
		//reset measured peak frequency values in chemicalShift.atoms
		for (i = 0; i< this.atoms.size(); i++) {
			this.atoms.get(i).resetMeasuredFrequencyValu();
			this.atoms.get(i).resetMeasFSpec();
		}
			
		for (specId = 0; specId < this.allSpectra.size(); specId++) {
			for (i = 0; i < this.allSpectra.get(specId).getAtomSize(); i++) {
				this.allSpectra.get(specId).getAtom(i).resetMeasuredFrequencyValu();
				this.allSpectra.get(specId).getAtom(i).resetMeasFSpec();
			}
		}
		
		
		
	}

	

	public void printMatrixOneByOne(double[][] matrix){
		int i, j, row, column;
		
		row = matrix.length;
		if (row == 0)
			System.out.println("The row size of matrix is zero, so the matrix cannot be printed.");
		else {
			column = matrix[0].length;
			if (column == 0)	
				System.out.println("The column size of matrix is zero, so the matrix cannot be printed.");
			else {
					for (i = 0; i< row; i++){
					for (j = 0; j < column; j++){
						System.out.println("matrix["+(i+1)+"]["+(j+1)+"]= "+ matrix[i][j] );
					}
				}
			}
		}
	}
	
	public void printMatrix(double[][] matrix){
		int i, j, row, column;
		
		row = matrix.length;
		if (row == 0)
			System.out.println("The row size of matrix is zero, so the matrix cannot be printed.");
		else {
			column = matrix[0].length;
			if (column == 0)	
				System.out.println("The column size of matrix is zero, so the matrix cannot be printed.");
			else {
				for (i = 0; i< row; i++){
					for (j = 0; j < column; j++){
						System.out.print(matrix[i][j] + " " );
					}
					System.out.println();
				}
			}
			System.out.println();
		}
	}
	/*
	 * This function returns the expected peak size of the whole experiment, by adding the expected peaks of each spectrum.
	 */
	private int getExpectedPeakSize() {
		int i;
		int e = 0;
		for (i = 0; i < this.allSpectra.size(); i++){
			e += this.allSpectra.get(i).getExpectedPeaksSize(); 
		}	
		return e;
	}


	/*
	 * This function is for extend spectrum. Constraint number is the atom number of the whole experiment (from all expected peaks)
	 * so, it is equal to the size of ChemicalShiftAssignment.atoms
	 */
	private void setObjectiveConstraintCounts(){
		int i;
		this.nObjectives = 4; //it can be increased directly here and also 2 other places: 1)final assignment in evaluateObjectives() to be updated, 2) in main function(2 main functions) negate for maximization.

		this.nConstraints = 0;// this.atoms.size(); //it can be increased directly here, and evaluateConstraints()  to be updated. And if 0:comment out evaluateConstraints call from evaluate method
		System.out.println("-number of objectives: "+ nObjectives);
		System.out.println("-number of constraints: "+ nConstraints);
		
		for (i= 0; i< this.allSpectra.size(); i++){
			this.allSpectra.get(i).setObjectiveConstraintCounts(this.nObjectives, this.nConstraints);
		}
		
	}
	

	/*
	 * This function creates a new solution according to the defined size and range of the expected peaks.   
	 */
	public Solution newSolution() {
		int s,i, j; 
		int mpSize;
		Solution solution = new Solution(nVariables, nObjectives, nConstraints); 
		
		i = 0;
		for (s = 0; s < allSpectra.size(); s++){
			mpSize = this.allSpectra.get(s).getMeasuredPeaksSize();
			for (j = 0; j < this.allSpectra.get(s).getExpectedPeaksSize(); j++){
				solution.setVariable(i,new BinaryIntegerVariable(-1, mpSize-1));
				i ++;
			}			
		}
		return solution;
	}
	
	
	/*
	 * evaluate method after extend to many spectra
	 * 
	 */
	public void evaluate(Solution solution){
		
		double G;
	 	ArrayList<Atom> atoms = new ArrayList<Atom>();
		int[] assignment = EncodingUtils.getInt(solution);
		//the order of the evaluateObjectives and evaluateConstraints is important since there are some methods in evaluateObjectives which are needed for the second one
	
		double[] objectives = new double[nObjectives]; 
		objectives = evaluateObjectives(assignment);
		G = this.currentG;
		atoms = this.atoms;

		solution.setObjectives(objectives);

		solution.setG(G);
		solution.setAtoms(atoms);
		
		
	}


	/*
	 * This function is to evaluate the score function values for the problems objectives.
	 * As input it receives the current assignment
	 * As output it returns the single objective values. 
	 */
	private double[] evaluateObjectives(int[] asg) {
	
		int i, j, expectedPeakSize, spectrumId, nextSpecStart, atomId, c;
		int expPeakIndex, measuredPeakIndex, expPeakIndexChemShift;
		int dimIndex = -1;
		double[] objectives = new double[nObjectives];
		double probabilityCondition, partTwo, partOne, G, a,b, d, x,y,z;
		double currentFrequency, deviation, devResult, difference, atomFrequency, thisW2, this_bonus, bonus, thisW1, thisMaxScore;
		double a_normalized, b_normalized, c_normalized, d_normalized;
		Peak expPeak = new Peak();
		Peak measuredPeak = new Peak();
	

		this_bonus = 0;
		a = 0; //shift normality-- to be maximized
		b = 0; //alignment-- to be maximized 
		d = 0; //degeneracy-- to be minimized 
		c = 0; //completeness -- to be maximized
		
	
		parseAssignment (asg);
		this.calculateDegeneracy(asg);
		
		thisW2 = 1/(-Math.log(cdf(-weightTwo,1)*2));
		thisW1 = 4;



		c = findAssignedAtomsOfExperiment(asg);
		for (spectrumId = 0; spectrumId < this.allSpectra.size(); spectrumId++){
			for (i = 0; i < this.allSpectra.get(spectrumId).getAtoms().size(); i++){
				this.allSpectra.get(spectrumId).getAtoms().get(i).clearMeasuredFrequency();
				this.allSpectra.get(spectrumId).calculateAverageFrequencyOfAtom(this.allSpectra.get(spectrumId).getAtoms().get(i), this.allSpectra.get(spectrumId).getAssignment(), i); 
			}
		}
		
		double score = 0;
	
		for (i = 0; i < this.atoms.size(); i++){ //for all atoms of the experiment

			deviation = 0;
			this.atoms.get(i).setDeviationValue(0);
			this.atoms.get(i).clearMeasuredFrequency();
			if (this.atoms.get(i).getAssigned(this.assignment.getAssignmentWhole()) == true){ //if this atom exists in an assigned expected peak
				Atom at = new Atom(this.atoms.get(i));
		
				expectedPeakSize = at.getExpectedPeakSize();
				partTwo = 0;
				for (j = 0; j < expectedPeakSize; j++){ //for all expected peaks of the atom
					expPeak = at.getExpectedPeak(j);
					partTwo = 0;
					spectrumId = expPeak.getSpectrumIndex();
			
					if (expPeak.getAssigned() == true){ //peak is assigned
						
						expPeakIndexChemShift = expPeak.getAssignmentIndex();
						measuredPeakIndex = asg[expPeakIndexChemShift];
						
						expPeakIndex = expPeak.getIndex();
						
						measuredPeak = this.allSpectra.get(spectrumId).getMeasuredPeaks().get(measuredPeakIndex);
						dimIndex = at.getExpectedPeaksIndexMapperOf(j);
						currentFrequency = measuredPeak.getChemicalShift(dimIndex);
						this.atoms.get(i).addMeasuredFrequency(currentFrequency);
						atomId = at.getIndicesSingleSpectrum(spectrumId); //index of that atom in Spectrum.atoms					
				
						deviation = Math.abs(currentFrequency - this.allSpectra.get(spectrumId).getAtoms().get(atomId).getAverageFrequency());
						devResult =  this.atoms.get(i).getAccuracy()/(this.maximumDeviation);
									
						probabilityCondition = cdf(-1*deviation, devResult)*2;
							
						probabilityCondition = Math.max(probabilityCondition, smallestValue);
						partTwo += ((Math.log(probabilityCondition)*this.weightThree) + 1)/this.allSpectra.get(spectrumId).getDegeneracyValue(expPeakIndex);						
						b += (Math.log(probabilityCondition)*this.weightThree) + 1;
						
						y = this.allSpectra.get(spectrumId).getDegeneracyValue(expPeakIndex); // these lines are used to get the d result below as double value without casting 
						x = 1/y;
						d += y;  //degeneracy scores that need to be maximized (it was d += x before; but changed now to y)
						
						this.atoms.get(i).increaseDeviationValue(partTwo);
						
				} else {					//Peak is not assigned

						
					}
				}

			
				//partFreq calculation starts here
       
				this.atoms.get(i).calculateAverageFrequency();

				if (this.atoms.get(i).getCsTableFrequency() > 0) { 
					difference = Math.abs(this.atoms.get(i).getAverageFrequency()- this.atoms.get(i).getReferenceFrequency());
					atomFrequency = cdf(-1.0*difference,this.atoms.get(i).getCsTableFrequency())*2;
					atomFrequency = Math.max(atomFrequency, smallestValue);
					
					partOne = Math.log(atomFrequency)*thisW2 + 1;
					a += partOne;
					
					if (this_bonus > 0 ){
						bonus = Math.exp(-0.5*(Math.pow(difference/this.atoms.get(i).getAccuracy(),2)))*this_bonus;  
						partOne = partOne + bonus;
					}
					
				} else{
					partOne = 0;
				}
			
				this.atoms.get(i).setFrequencyValue(partOne*thisW1);
				this.atoms.get(i).setGlobalScore (this.atoms.get(i).getDeviationValue() + this.atoms.get(i).getFrequencyValue());
				score += this.atoms.get(i).getGlobalScore();
							
				
			} else {  //atom is not assigned
				this.atoms.get(i).setAverageFrequency(0); //average frequency of that atom is 0.
			}
			
		}


		G = score / this.maximumScore;
		this.currentG = G;

		objectives [0] = -a;
		objectives [1] = -b;
		objectives [2] = d;
		objectives [3] = -c;
	
		return objectives;
	}

	/*
	 * This function sets the count of the individuals of the algorithm's population as given input value. 
	 */
	public void setPopulationSize(int value) {
		this.populationSize = value;
	}
	public int getPopulationSize() {
		return this.populationSize;
	}

	/*
	 * This method fills the constraints[] and returns it
	 */
	public double[] evaluateConstraints(int[] asg){
		double[] constraints = new double[nConstraints];
		int i, j ,mfSize, a,b;
		double accuracy, difference;
		double punishment = 0;
		
		
		/*
		 * Evaluates the value of the first constraint, i.e. variations of chemical shifts assigned to an atom.
		 */
		for (i = 0; i < this.atoms.size(); i++){ // for all atoms of the spectrum
			punishment = 0;
			mfSize = this.atoms.get(i).getMeasuredFrequencySize();
			accuracy =  this.atoms.get(i).getAccuracy();
			if (mfSize > 0){  // if there is at least one measured frequency assigned to this atom
			
				for (a = 0; a < (mfSize-1); a++){
					for (b = a+1; b < mfSize; b++){
						difference = Math.abs(this.atoms.get(i).getMeasuredFrequency(a)- this.atoms.get(i).getMeasuredFrequency(b));
						if (difference > accuracy){ // add punishment since the difference exceeds the given tolerance,i.e. accuracy
							punishment += difference;
						}
					}
				}
			}
			constraints[i] = punishment;
		}
		
		return constraints;
	}
	
	
	/*
	 * This function is for extension to many spectrum
	 */
	private int findAssignedAtomsOfExperiment(int[] assignment) {
		int i, j, k, index, atomsIndex, spectrumId, expectedPeakSize, dimension, expPeakIndex, startIndex, assignmentIndex, counter;
		String atomName;
		resetAtomAssignments();
		counter = 0;

		
		index = 0; 
		
		startIndex = 0;
		for (spectrumId = 0; spectrumId < this.allSpectra.size(); spectrumId++){
			expectedPeakSize = this.allSpectra.get(spectrumId).getExpectedPeaksSize(); //asg size of that spectrum
			dimension = this.allSpectra.get(spectrumId).getDimension();
			Atom[] assignedAtoms = new Atom[dimension];
			
			for (i = 0; i < dimension; i++){
				assignedAtoms[i] = new Atom();
			}
			for (j = 0; j < expectedPeakSize; j++){
				expPeakIndex = this.allSpectra.get(spectrumId).getExpectedPeaks().get(j).getIndex();  //shows the index of this expected peak in Spectrum.expectedPeaks; so it is the index of that exp peak in the assignemtn of the spectrum
				
				if (assignment[index] != -1){ //if index. expected peak is assigned
					counter ++;
					this.allSpectra.get(spectrumId).getExpectedPeaks().get(j).markAssigned();
					assignedAtoms = this.allSpectra.get(spectrumId).getExpectedPeaks().get(j).getAtoms();
					for (k = 0; k < dimension; k++){
						atomName = assignedAtoms[k].getName();
						atomsIndex = assignedAtoms[k].getAtomsIndex();
						this.atoms.get(atomsIndex).setAssigned(true);
					}

				} else if (assignment[index] == -1){ //peak is not assigned

				}
				index++;
			}
			
			startIndex += this.allSpectra.get(spectrumId).getExpectedPeaksSize();
			
		}
		
		for (i = 0; i < this.atoms.size(); i++){  //for all atoms in ChemicalShiftAssignment.atoms
			for (j = 0; j < this.atoms.get(i).getExpectedPeakSize(); j++ ){		//for all expected peaks of this atom ChemicalShiftAssignment.Atoms.expectedPeaks
				assignmentIndex = this.atoms.get(i).getExpectedPeak(j).getAssignmentIndex();
				
				if (assignment[assignmentIndex] ==-1){  //not assigned
					this.atoms.get(i).getExpectedPeak(j).setAssigned(false);
					for (k = 0; k < this.atoms.get(i).getExpectedPeak(j).getDimension(); k++)
						this.atoms.get(i).getExpectedPeak(j).getAtom(k).setAssigned(false);
					
				} else if (assignment[assignmentIndex] >= 0){  //assigned
					this.atoms.get(i).getExpectedPeak(j).setAssigned(true);
					for (k = 0; k < this.atoms.get(i).getExpectedPeak(j).getDimension(); k++)
						this.atoms.get(i).getExpectedPeak(j).getAtom(k).setAssigned(true);
				}
			}
		}

		return counter;
	}

	/*
	 * This function is to reset all atom assignments. 
	 */
	public void resetAtomAssignments(){
		int i,j;
		int atomCount = this.atoms.size();
		int epCount;
		
		for (i = 0; i < atomCount; i++){
		
			this.atoms.get(i).setAssigned(false);
			this.atoms.get(i).clearMeasuredFrequency();
			epCount = this.atoms.get(i).getExpectedPeakSize();
			for (j = 0; j < epCount; j++){
				this.atoms.get(i).getExpectedPeak(j).setAssigned(false);
			}
		}
		
		for (i = 0; i < this.allSpectra.size(); i++) {
			this.allSpectra.get(i).resetAtomAssignments();
		}
	}

	public double cdf(double x, double d){
		double cdf, e;
		Derf a = new Derf();
				
		e = x/(Math.sqrt(2.0)*d);
		cdf = 1.0 - 0.5*(a.erfcc(e));
		
		return cdf;
	}
	

	
	 /*
	  * This function is written to calculate degeneracy for extension to many spectrum 
	  */
	private void calculateDegeneracy(int[] asg) {

		int i;
		for (i = 0; i < this.allSpectra.size(); i++){
			this.allSpectra.get(i).calculateDegeneracy(this.allSpectra.get(i).getAssignment());
		}
	}

	/*
	 * This function is for extension to many spectrum. It parses the whole assignment and creates its own assignment in each spectrum
	 */
	private void parseAssignment(int[] wholeAssignment) {
		int i, firstIndex,asgSize,specIndex;
		

		for ( specIndex = 0; specIndex < this.allSpectra.size(); specIndex++){
			asgSize = this.allSpectra.get(specIndex).getAssignmentSize();
			int[] asg = new int [asgSize];
			firstIndex = 0;
			for (i = 0; i < specIndex; i++){
				firstIndex += this.allSpectra.get(i).getAssignmentSize();
			}
			for (i = 0; i < this.allSpectra.get(specIndex).getAssignmentSize(); i++){
				asg[i] = wholeAssignment[firstIndex + i];
			}
			
			this.allSpectra.get(specIndex).setAssignment(asg);
			
			
		}
	}

	/*
	 * This function initializes the problem model including the whole library, reference frequency values, atoms, expected and measured peaks and spectra. 
	 */
	private void initialize() throws IOException {
		
		int i; 
		
		// creates library in all spectra
		for (i = 0; i < this.allSpectra.size(); i++){
			this.allSpectra.get(i).setLibrary(this.library);
			this.allSpectra.get(i).setAminoAcidSequence(this.aminoacidSequence);
		}
		
		//assigns atom reference frequency values in all spectra 
		for (i = 0; i < this.allSpectra.size(); i++){
			this.allSpectra.get(i).assignReferenceFrequency();
		}
	
		createAtoms();
		setInitialValues();
		
		findExpectedPeakIdsInEachSpectrum();
		
		initData();
		
	}
	
	/*
	 * This function triggers the initialization process for the score values for the initialization. 
	 */
	private void initData(){
		
		initScore();
		
	}

	
	/*
	 * This function compares the given two arrays and returns true if they are the same,
	 * it returns false, otherwise.
	 */
	private boolean compareTwoBooleanArrays(boolean[] array1, boolean[] array2) {
		int i;
		if (array1.length!= array2.length)
			return false;
		
		for (i = 0; i < array1.length; i++) {
			if (array1[i] != array2[i])
				return false;
		}
		return true;
	} 
	
	/*
	 * This function sets the given two expected peaks as equivalent peaks for each other. 
	 */
	private void addEquiv(int specId, int peak1, int peak2) {
		
		this.allSpectra.get(specId).getExpectedPeak(peak1).addToEqualPeaks(peak2);
		this.allSpectra.get(specId).getExpectedPeak(peak2).addToEqualPeaks(peak1);
	
	}

	/*
	 * This function initializes the score related variables. 
	 */
	private void initScore(){
		
		int specCount = this.allSpectra.size();
		double[] specPExist = new double[specCount];

		
		specPExist = this.countProbs();
		this.setPExistSum(specPExist);
		
		this.countAtomsPerDim();
		this.calculateAtomsValues();
		calculateEliteAssignment();

	}

	/*
	 * This function saves the elite assignment. 
	 */
	private void calculateEliteAssignment() {
		
		int atomId, atomsIndexInSpec, specId, peakId, expIdInSpec, idim;
		double w1 = 0;
		double maxScore = 0;
		w1 = 0;	
		this.thisMaximumScore = 0;
		for (atomId = 0; atomId < this.atoms.size(); atomId++) {
			if (this.atoms.get(atomId).getExpectedPeakSize()>0) {
				this.atoms.get(atomId).setAverageScoreLocal(0.5);
				this.atoms.get(atomId).setMeanDev(0.2); 
				this.thisMaximumScore = this.thisMaximumScore + this.thisWeightOne;
				w1 = w1 + this.thisWeightOne;
				peakId = getNextPeak(atomId, true);
				
				while (peakId >= 0 ) {
					specId = this.expectedPeakSpecId[peakId];
					expIdInSpec = this.expectedPeakExpIdInSpec[peakId];
					
					for (idim = 0 ; idim < this.allSpectra.get(specId).getDimension(); idim ++) {
						if (this.allSpectra.get(specId).getExpectedPeak(expIdInSpec).getAtom(idim).getName().equals(this.atoms.get(atomId).getName())) {
							this.thisMaximumScore = this.thisMaximumScore  + 1;
						}
					}
					peakId = getNextPeak(atomId, false);
					
				}
			    
			}
		}
	}
	
	/*
	 * This function calculates the number of different atoms for each spectrum's dimension.
	 */
	private void countAtomsPerDim(){
		
		int ispec, idim, atomCounter, atomId, peakId;
		boolean foundP = false;
		int nrSpectra, specId, expId;
		atomCounter = 0;
		
		nrSpectra = this.allSpectra.size();

		for (ispec = 0; ispec < nrSpectra; ispec++){
			this.allSpectra.get(ispec).getMinMaxFreq();
			for (idim = 0; idim < this.allSpectra.get(ispec).getDimension(); idim++){
				atomCounter = 0;
				for (atomId = 0; atomId < this.atoms.size(); atomId++){
					foundP = false;
					peakId = getNextPeak(atomId, false);
					while (peakId != -1){
						expId = this.expectedPeakExpIdInSpec[peakId];
						specId = this.expectedPeakSpecId[peakId];
						if (ispec == specId){
							
							if (atomId == this.allSpectra.get(specId).getExpectedPeak(expId).getAtom(idim).getAtomsIndex()){
								foundP = true;
							}
						}
						peakId = getNextPeak(atomId, false);
	             	}
	                if (foundP){
		                 atomCounter++;
	                }
				}
				this.allSpectra.get(ispec).setNumberOfAtoms(idim, atomCounter);
			}

		}
	}
	
	/*
	 * This function sets the pExistSum, degeneray factor and pValues for each spectrum.
	 */
	private void setPExistSum(double[] vals){
		int ispec;
	    int specCount = this.allSpectra.size();
		
	    for (ispec = 0; ispec < specCount; ispec++){
	    	this.allSpectra.get(ispec).setpExistSum(vals[ispec]);
	    	this.allSpectra.get(ispec).setDegeneracyFactor((double)this.allSpectra.get(ispec).getMeasuredPeaksSize()/vals[ispec]);
	    	this.allSpectra.get(ispec).setPValue(Math.min(vals[ispec]/(double)this.allSpectra.get(ispec).getMeasuredPeaksSize(), 1.0));
	    }
	    
	}
	
	/*
	 * This function calculates the sum of all probabilities for each spectrum. 
	 */
	private double[] countProbs(){
		int specCount = this.allSpectra.size();
		int i,spec;
		double total = 0;
		double pExist; //probability of the expected peak
		double[] specProbs = new double[specCount];
		
		for (i = 0; i < specCount; i++){
			specProbs[i] = 0;
		}
		
		for (spec = 0; spec < specCount; spec++){
			for (i = 0; i < this.allSpectra.get(spec).getExpectedPeaksSize(); i++){
				pExist = this.allSpectra.get(spec).getExpectedPeak(i).getPeakProbability();
				specProbs[spec] += pExist;
				total += pExist;
			}
		}
		
		return specProbs;
	}

	/*
	 * This function sets the scaling factor values for each atom in each spectrum.
	 */
	private void calculateAtomsValues(){
	     int ispec, type, minPeaksInSpec, pSpecType = 0, atomId, atomCount, i, peakId;
	     double[] sumTypePeaks;
	     double[] relSpecType;
	     double[] weight;
	     double minExpPeaks, sumWeight, sumFacExist, pExists;
	     int spectrumCount = this.allSpectra.size(), expId, specId;
	     
	     sumTypePeaks = new double[this.getNrSpecTypes()];
	     relSpecType = new double[this.getNrSpecTypes()];
	     weight = new double[this.getNrSpecTypes()];
	     
	     
	     for (ispec = 0; ispec < this.getNrSpecTypes(); ispec++){
	    	 sumTypePeaks[ispec] = 0;
	    	 relSpecType[ispec] = 0;
	     }
	 
	     for (ispec = 0; ispec < spectrumCount; ispec++){
	    	 	type = this.allSpectra.get(ispec).getCategory();
	    	 	sumTypePeaks[type] += this.allSpectra.get(ispec).getExpectedPeakCount();
	    	 	if (ispec == 0){
	    	 		minExpPeaks = this.allSpectra.get(ispec).getExpectedPeakCount();
	    	 		minPeaksInSpec = this.allSpectra.get(ispec).getExpectedPeaksSize();
	    	 	}
	     }
	     
	     atomCount = this.atoms.size();
	     for (atomId = 0; atomId < atomCount; atomId++){
	    	 this.atoms.get(atomId).initializeScalingFactor(this.allSpectra.size(), this.getNrSpecTypes());
	    	 sumWeight = 0;
	    	 for (i = 0 ; i < this.spectrumWeight.length; i++){
	    		 weight[i] = this.spectrumWeight[i];
	    		 sumWeight += weight[i];
	    	 }
	    	 if (this.weightTypeValue != 0){
	    		 
	    	 } else {
	    		 for (i = 0; i < relSpecType.length; i++){
	    			 relSpecType[i] = weight[i] / sumWeight / sumTypePeaks[i];
	    		 }
	    	 }
	    	 sumFacExist = 0;
	    	 peakId = getNextPeak(atomId, true);
	    	 while (peakId != -1){
	    		 expId = this.expectedPeakExpIdInSpec[peakId];
	    		 specId = this.expectedPeakSpecId[peakId];
	    		 pExists = this.allSpectra.get(specId).getExpectedPeak(expId).getpExists();
	    		 this.atoms.get(atomId).increaseProbabilitySum(pExists);
	    		 this.atoms.get(atomId).increaseExpInSpec(specId, pExists);
	    		 sumFacExist += pExists * relSpecType[this.allSpectra.get(specId).getCategory()];
	    		 peakId = getNextPeak(atomId, false);
	    	 }
	    	 this.atoms.get(atomId).setScalingFactorSpectrum(relSpecType, sumFacExist);
		}
	}
	
	/*
	 * This function sets the expected peak Ids in each spectrum.
	 */
	private void findExpectedPeakIdsInEachSpectrum() {
		int i, expIdSpect, specId;
		int results[] = new int[2];
		this.expectedPeakExpIdInSpec = new int[this.getExpectedPeakSize()];
		this.expectedPeakSpecId = new int[this.getExpectedPeakSize()];
		
		expIdSpect = -1;
		specId = -1;
		
		for (i = 0; i < this.getExpectedPeakSize(); i++){
			results = findExpIdInSpectrum(i);
			specId = results[0];
			expIdSpect = results[1];
			this.expectedPeakSpecId[i] = specId;
			this.expectedPeakExpIdInSpec[i] = expIdSpect;
		}
		
	}


	
	/*
	 * This function receives the expectedPeakId from ChemicalShiftAssignment.assignment and returns the expected Peak Id of that expected peak in its own spectrum
	 * results[0] = specId; results[1]= expId
	 */
	private int[] findExpIdInSpectrum(int expId) {
		int firstIndex,lastIndex, asgSize,specId, expIdSpect;
		int results[] = new int[2];
		expIdSpect = -1;
		specId = -1;


		firstIndex = 0;
		for (specId = 0; specId < this.allSpectra.size(); specId++){
			asgSize = this.allSpectra.get(specId).getAssignmentSize();
			lastIndex = firstIndex + asgSize - 1;  
			if (expId <= lastIndex){
				expIdSpect = expId - firstIndex;
				break;
			}

			firstIndex += this.allSpectra.get(specId).getAssignmentSize();
		}
		results[0] = specId;
		results[1] = expIdSpect;
		return results;
	}

	
	/*
	 * This function creates the atoms of the whole experiment and assigns the expected peaks to them. 
	 * All atoms have been already created in each spectrum, this method scans these atoms and creates the atoms in chemicalshiftassignment class.
	 * This is static for the experiment. Assignment affects are not evaluated yet.
	 *
	 */
	public void createAtoms() throws IOException{

		int i, j, atomIndex, spec, a, startIndex, residueNumber;
		String atomName, genericAtomName; 
		Peak pk;
		Atom at; 

		//atoms of the chemical shift assignment experiment, so all spectrum's atoms 
		try{
			startIndex = 0; //index of the first value of this spectrum in the chemicalshiftassignment.assignemnt
			for (spec = 0; spec < this.allSpectra.size(); spec++){ //for all spectrum
				for (a = 0; a < this.allSpectra.get(spec).getAtomSize(); a++){ //for all atoms of that spectrum
					atomName = this.allSpectra.get(spec).getAtom(a).getName();
					atomIndex = atomExists(atomName);
					if (atomIndex >= 0){ // atom exists; so add its values to the existing atom
						this.atoms.get(atomIndex).setIndicesSingleSpectrum(spec, a);
						
					} else if (atomIndex == -1){ //atom is new; so create new atom, add its values and add to list
						residueNumber = Integer.parseInt(atomName.substring(atomName.indexOf(".")+1));
						genericAtomName = atomName.substring(0, atomName.indexOf("."));
						at = new Atom(atomName, this.allSpectra.size(), spec, a, residueNumber, genericAtomName);
						at.setIndex(this.atoms.size());
						at.setReferenceFrequency(this.allSpectra.get(spec).getAtom(a).getReferenceFrequency());
						at.setCSTableFreq(this.allSpectra.get(spec).getAtom(a).getCsTableFrequency());
						atomIndex = this.atoms.size(); //if there were 2 atoms, after adding this new atom, its index will be 2 since list starts from 0.
						this.atoms.add(at);
						
					}  else{
						throw new IOException("Atom index is not found for expected peak in createAtoms method of Spectrum.class");
					}	
					
					//add all expected peaks of this atom from spectrum to ChemicalShiftAssignment.atoms.expectedPeaks
					this.allSpectra.get(spec).copyExpectedPeaksFromOneAtomToAnotherAtom(this.allSpectra.get(spec).getAtom(a), this.atoms.get(atomIndex), startIndex);
				}
			startIndex += this.allSpectra.get(spec).getExpectedPeaksSize();
			}
			setAtomIndicesInExpectedPeaks();

		} finally {
	
		}
	

	}



	/*
	 * This function is designed for extend to many spectrum. It sets the index of the atom in the ChemicalShiftAssignment.atoms in the Spectrum.expectedPeaks
	 * 
	 */
	private void setAtomIndicesInExpectedPeaks() {
		int spec, i, j, dimension, atomsIndex;
		String atomName;
		Peak pk = new Peak();
		
		for (spec = 0; spec < this.allSpectra.size(); spec++){ //for all spectra
			dimension = this.allSpectra.get(spec).getDimension();
			for (i = 0; i < this.allSpectra.get(spec).getExpectedPeaksSize(); i++) {
				for (j = 0; j < dimension; j++){
					pk = this.allSpectra.get(spec).getExpectedPeaks().get(i);
					atomName = pk.getAtom(j).getName();
					atomsIndex = findAtomsIndex (atomName);
					if (atomsIndex == -1)  //The atom does not exist 
						System.out.println("ChemicalShiftAssignment.findAtomsIndex: Atom of spectrum does not exist in atoms array but it is supposed to be.");
					else if (atomsIndex >= 0){ //atom exists
						this.allSpectra.get(spec).getExpectedPeaks().get(i).getAtom(j).setAtomsIndex(atomsIndex);
					}
				}
				
				
				}
			}	
	}


	/*
	 * This function is designed for extend spectrum. Finds the ChemicalShiftAssignment.atoms index of the given atom
	 */
	private int findAtomsIndex(String atomName) {
		int atomsIndex = -1;
		int i;
		
		for (i = 0; i < this.atoms.size(); i++){
			if (atomName.equals(this.atoms.get(i).getName())){
				atomsIndex = i;
				return atomsIndex;
			}
		}
		
		return atomsIndex;
	}

	/*
	 * This function search whether the input atom already exist in the experiment.
	 * It returns -1 if atom does not exist in the atoms 
	 * returns the index of the atom otherwise
	 * */
	int atomExists(String atomName){
		int index = -1;
		int i;
		
		for (i = 0; i < this.atoms.size(); i++){
			if (atomName.equals(this.atoms.get(i).getName())){
				index = i;
				return index;
			}
		}
		return index;
	}

	/*
	 * This function sets the initial values for spectrum variables.
	 */

	private void setInitialValues() {
		setAccuracyValues();  
		maximumDeviation = 4;
		initializeMaxScore();
		this.spectrumWeight = new double[1];
		this.spectrumWeight[0] = 1;
	}

	/*
	 * This function sets the accuracy values for each atom of the NMR experiment.
	 */
	public void setAccuracyValues(){
		int i, size;
		size = this.atoms.size();
		
		for (i = 0; i < size; i++){
			this.atoms.get(i).setAccuracyValues();
		}
	}
	
	/*
	 * This function calculates the denominator for the Global score calculation that will be used for normalization.
	 */
	private void initializeMaxScore() {
		
		int i;
		double denom = 0;
		double w1_a = 4; 
		double w2_an = 1; 
		double n_atoms;
		for (i = 0; i < this.atoms.size(); i++){
			n_atoms = this.atoms.get(i).getExpectedPeakSize();
			denom += (w1_a +(w2_an*n_atoms));
		}
		this.maximumScore = denom;
	}

	/*
	 * This function loads input files: sequence fiel, shift reference file and the library file.
	 */
	private void loadFiles()throws IOException{
		int i;
		String spectrum;
		File file;
		try{
			for(i = 0; i< this.spectrumNames.length; i++){
				spectrum = path + spectrumNames[i]; //path+HNCA
				spectrum= spectrum.replaceAll("\\s","");
				loadSpectrum(spectrum,i);
				System.out.println(spectrumNames[i]+" spectrum ("+this.allSpectra.get(i).getMeasuredPeaksSize()+" measured, "+this.allSpectra.get(i).getExpectedPeaksSize()+" expected peaks) is read. ");
			}
			loadSequenceFile();
			System.out.println("Sequence file '"+sequenceFile+"' ("+ aminoacidSequence.size() +" aminoacids) is read.");
			
			
			loadShiftReferenceFile();
			System.out.println("Shift reference file (real values) '"+shiftReferenceFile+"' ("+this.shiftReferences.size()+" atoms) is read.");
			
			loadLibraryFile();
			System.out.println("Library file '" +libraryFile+ "' ("+this.library.getAtomTypeCount()+" atom types and "+this.library.getResidueCount()+" residues) is read.");
			
		} finally {
			
		}
	}
	
	/*
	 * This function prints the given solution candidate.
	 */
	private void printIndividual(Solution individual) {
		int i;
		
		for (i = 0; i < individual.getNumberOfVariables(); i++) {
			System.out.print(individual.getVariable(i) + " ");
		}
		System.out.println();
		
	}

	/*
	 * This function returns a random variable.
	 */
	private double getRandomNumber() {
		double rand;
		
		rand = Math.random();
		
		return rand;		
	}
	
	/*
	 * This function loads the library file and initializes all related attributes and objects. 
	 */
	private void loadLibraryFile() throws IOException{
		String fileName = path + libraryFile;
		int atomTypeCount,i, j, atomTypeNumber, hydrogenCode, order;
		AtomType at;
		Residue re;
		CSTableReference csReference;
		String atomType;
		float coreRadius, lastValue;
		String name;
		int angleCount, atomCount, firstAtomNumber, lastAtomNumber;
		int dihedralNumber,dihedralLastAtom;
		String dihedralName;
		int[] dihedralAtoms = new int[4]; // the number of 4 atoms that define the dihedral angle
		ArrayList<DihedralAngle> dihedralAngles = new ArrayList<DihedralAngle>();
		DihedralAngle da = new DihedralAngle();
		ArrayList<Atom> atoms = new ArrayList<Atom>();
		Atom atom = new Atom();

		CSTableReference cs = new CSTableReference();
		
		int atomNumber, atomPseudoAtomNumber, cSTableCount;
		String atomName, atom_type;
		float[] atomCoordinates = new float[3]; // x-, y- and z-coordinates 
		int[] atomCovalentConnectivityAtoms = new int[4]; //	four atom numbers indicating covalent connectivities (if there are less than four connectivities, the corresponding numbers are set to 0)
		int idNumber, notUsedNumber;
		String residueName, cSAtomName;
		double referenceFrequency, second, third, fourth;
		
		Pattern atomTypesHeader = Pattern.compile("ATOMTYPES( |\t)+(\\d+)(.*?)");
		Pattern atomTypeLine = Pattern.compile("( |\t)+(\\d+)( |\t)+(\\S+)*( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+)( |\t)+(\\d+)( |\t)+(-?\\d+\\.\\d+)(.*?)");
		Pattern residueLine = Pattern.compile("RESIDUE( |\t)+(\\S+)*( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)(.*?)");
		Pattern angleLine = Pattern.compile("( |\t)+(\\d+)( |\t)+(\\S+)*( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(.*?)");
		Pattern atomLine = Pattern.compile("( |\t)+(\\d+)( |\t)+(\\S+)*( |\t)+(\\S+)*( |\t)+(\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)"
				+ "( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)( |\t)+(\\d+)(.*?)");
		Pattern csTableHeader = Pattern.compile("CSTABLE( |\t)+(\\d+)(.*?)");
		Pattern csTableLine = Pattern.compile("( |\t)+(\\d+)( |\t)+(\\S+)*( |\t)+(\\S+)*( |\t)+(\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)(.*?)");
		 
		CommentedLineReader lineReader = null;
		String line = null;
		Matcher matcher = null;
		Matcher matcherCSTable = null;
		
		try{
			File file = new File(fileName);
			if (file.exists() == false){
				throw new IOException("Library file "+fileName+" does not exist!");
			}
			Reader reader = new FileReader(file);
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); //the number of dimension line
			
			matcher = atomTypesHeader.matcher(line);
			if (matcher.matches()) {
				atomTypeCount = Integer.parseInt(matcher.group(2));
				this.library.setAtomTypeCount(atomTypeCount);
							
			} else {
				throw new IOException("Library file not properly formatted: invalid 'ATOMTYPES header' line");
			}
			
			for (i = 0; i < atomTypeCount; i++){
				line = lineReader.readLine(); //the number of dimension line
				matcher = atomTypeLine.matcher(line);
				if (matcher.matches()) {
					atomTypeNumber = Integer.parseInt(matcher.group(2));
					atomType = matcher.group(4);
					coreRadius = Float.parseFloat(matcher.group(6));
					hydrogenCode = Integer.parseInt(matcher.group(8));
					order = Integer.parseInt((matcher.group(10)));
					lastValue = Float.parseFloat((matcher.group(12)));
					
					at = new AtomType(atomTypeNumber, atomType, coreRadius, hydrogenCode, order, lastValue);
					this.library.addAtomType(at);
					
				} else {
					throw new IOException("Library file not properly formatted: invalid 'atom type' line");
				}
			}
			
			line = lineReader.readLine();
		
			while(line!=null){
				
				matcher = residueLine.matcher(line);
				
				while (matcher.matches() == false){
					line = lineReader.readLine();
					matcherCSTable = csTableHeader.matcher(line);
					if (line == null)
						 break;
					if (matcherCSTable.matches() == true){
						break;
					}
						
					matcher = residueLine.matcher(line);
				}
				
				if (matcher.matches()) {
					name = matcher.group(2);
					angleCount = Integer.parseInt(matcher.group(4));
					atomCount = Integer.parseInt(matcher.group(6));
					firstAtomNumber = Integer.parseInt(matcher.group(8));
					lastAtomNumber = Integer.parseInt(matcher.group(10));
					dihedralAngles.clear();
					for (i = 0; i<angleCount; i++){
						line = lineReader.readLine();
						matcher = angleLine.matcher(line);
						if (matcher.matches()) {
							
							dihedralNumber = Integer.parseInt(matcher.group(2));
							dihedralName = 	matcher.group(4);
							dihedralAtoms[0] = Integer.parseInt(matcher.group(12));
							dihedralAtoms[1] = Integer.parseInt(matcher.group(14));
							dihedralAtoms[2] = Integer.parseInt(matcher.group(16));
							dihedralAtoms[3] = Integer.parseInt(matcher.group(18));
							dihedralLastAtom = Integer.parseInt(matcher.group(20));
							
							da = new DihedralAngle(dihedralNumber, dihedralName, dihedralAtoms, dihedralLastAtom);
									
							dihedralAngles.add(da);
	
						} else {
							throw new IOException("Library file not properly formatted: invalid 'dihedral angle' line");
						}
					}
					atoms.clear();
					for (i = 0; i<atomCount; i++){
						line = lineReader.readLine();
						matcher = atomLine.matcher(line);
						if (matcher.matches()) {
							atomNumber = Integer.parseInt(matcher.group(2));
							atomName = matcher.group(4);
							atom_type = matcher.group(6);
							atomCoordinates[0] = Float.parseFloat(matcher.group(12));
							atomCoordinates[1] = Float.parseFloat(matcher.group(14));
							atomCoordinates[2] = Float.parseFloat(matcher.group(16));
							atomCovalentConnectivityAtoms[0] = Integer.parseInt(matcher.group(18));
							atomCovalentConnectivityAtoms[1] = Integer.parseInt(matcher.group(20));
							atomCovalentConnectivityAtoms[2] = Integer.parseInt(matcher.group(22));
							atomCovalentConnectivityAtoms[3] = Integer.parseInt(matcher.group(24));
							atomPseudoAtomNumber = Integer.parseInt(matcher.group(26));
					

							atom = new Atom (atomNumber, atomName, atom_type, atomCoordinates, atomCovalentConnectivityAtoms, atomPseudoAtomNumber);
							atoms.add(atom);

							
						} else {
							throw new IOException("Library file not properly formatted: invalid 'dihedral angle' line");
						}
					}
	
					re = new Residue(name, angleCount, atomCount, firstAtomNumber, lastAtomNumber, dihedralAngles, atoms);
					this.library.addResidue(re);
								
				}
				// start reading the CSTABLE part
				else if (matcherCSTable.matches() == true) {
				
					cSTableCount = Integer.parseInt(matcherCSTable.group(2));
			
					this.library.setCSTableCount(cSTableCount);
					line = lineReader.readLine();
					
					while(line!=null){
						
						matcherCSTable = csTableLine.matcher(line);
						
						if (matcherCSTable.matches() == true){
							idNumber = Integer.parseInt(matcherCSTable.group(2));
							residueName = matcherCSTable.group(4);
							cSAtomName = matcherCSTable.group(6);
							notUsedNumber = Integer.parseInt(matcherCSTable.group(8));
							referenceFrequency = Double.parseDouble(matcherCSTable.group(10));
							second = Double.parseDouble(matcherCSTable.group(12));
							third = Double.parseDouble(matcherCSTable.group(14));
							fourth = Double.parseDouble(matcherCSTable.group(16));
						
					
							
							cs = new CSTableReference(idNumber, residueName, cSAtomName, notUsedNumber, referenceFrequency, second, third, fourth);
							this.library.addReferenceToCSTable(cs);
							
						}
						line = lineReader.readLine();	
					}
					System.out.println("CS Table ("+ this.library.getcSTableCount()+" entries) is read. " );
					
			
				}
				
				else if (line!=null) {
					throw new IOException("Library file not properly formatted: invalid 'residue header' line");
				}
				
				line = lineReader.readLine();
				
			
			
			}
		} finally {
			if (lineReader != null) {
				lineReader.close();
			}
		}
	
		
		this.library.assignStatistics(this.aminoacidSequence);
		lineReader.close();
	}
	
	/*
	 * This function loads the shift reference file
	 */
	private void loadShiftReferenceFile() throws IOException{
		int i, atomNumber, residueNumber;
		float chemicalShift, chemicalShiftError;
		String atomName;
		String fileName = path + shiftReferenceFile;
		CommentedLineReader lineReader = null;
		String line = null;
		Matcher matcher = null;
		ShiftReference sr;
		
		Pattern shiftReferenceLine = Pattern.compile("( |\t)+(\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(-?\\d+\\.\\d+)( |\t)+(\\S+)*( |\t)+(\\d+)(.*?)");
		try{
			File file = new File(fileName);
			if (file.exists() == false){
				throw new IOException("Input peak file "+fileName+" does not exist!");
			}
			Reader reader = new FileReader(file);
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); //the number of dimension line
			
			while (line!=null){
				matcher = shiftReferenceLine.matcher(line);
				if (matcher.matches()) {
					
					atomNumber = Integer.parseInt(matcher.group(2));
					chemicalShift = Float.parseFloat(matcher.group(4));
					chemicalShiftError = Float.parseFloat(matcher.group(6));
					atomName = matcher.group(8);
					residueNumber = Integer.parseInt(matcher.group(10));
									
					sr= new ShiftReference(atomNumber, chemicalShift, chemicalShiftError, atomName, residueNumber);
					this.addShiftReference(sr);
										
				} else {
					throw new IOException("Shift reference file not properly formatted: invalid 'shift reference' line");
				}
				
				line = lineReader.readLine();			
			}
			
		} finally {
			if (lineReader != null) {
				lineReader.close();
			}
		}
		lineReader.close();
	}
	
	public void addShiftReference(ShiftReference sr){
		ShiftReference mySr = new ShiftReference(sr);
		this.shiftReferences.add(mySr);
	}
	
	/*
	 * This function loads the sequence file that defines the sequence of the amino acids in the NMR experiment. 
	 */
	private void loadSequenceFile() throws IOException{
		String fileName = path + sequenceFile;

		CommentedLineReader lineReader = null;
		String line = null;
		

		try{
			File file = new File(fileName);
			if (file.exists() == false){
				throw new IOException("Input peak file "+fileName+" does not exist!");
			}
			Reader reader = new FileReader(file);
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); //the number of dimension line
			
			while (line!=null){
				line = line.replaceAll("\\s+","");
				aminoacidSequence.add(line);
				line = lineReader.readLine();
			}
			
		} finally {
			if (lineReader != null) {
				lineReader.close();
			}
		}
		lineReader.close();
	}
	
	/*
	 * This function adds the given spectrum to the current problem.
	 */
	private void loadSpectrum(String spectrum, int index) throws IOException {

			
		this.allSpectra.add(new Spectrum(spectrum, index));
	}
	

	/*
	 * This function reads the names of the input files from the configuration file and loads the name of the input files accordingly.
	 */
	
	private void loadFileNames(Reader reader) throws IOException {
		int i;
		Pattern peakLine = Pattern.compile("peaks:=(.*?)");
		Pattern sequenceLine = Pattern.compile("sequence:=(.*?)");
		Pattern shiftReferenceLine = Pattern.compile("shiftreference:=(.*?)");	
		Pattern libraryLine = Pattern.compile("library:=(.*?)");
		
		CommentedLineReader lineReader = null;
		String line = null;
		Matcher matcher = null;
		
		try{
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); 
			matcher = peakLine.matcher(line);
			if (matcher.matches()) {
				this.spectrumNames = matcher.group(1).split(",");
			} else {
				throw new IOException("The configuration file not properly formatted: invalid 'peaks=' line");
			}
			
			line = lineReader.readLine(); 
			matcher = sequenceLine.matcher(line);
			if (matcher.matches()) {
				this.sequenceFile = matcher.group(1);	
			} else {
				throw new IOException("The configuration file not properly formatted: invalid 'sequence=' line");
			}
			
			line = lineReader.readLine(); 
			matcher = shiftReferenceLine.matcher(line);
			if (matcher.matches()) {
				this.shiftReferenceFile = matcher.group(1);
			} else {
				throw new IOException("The configuration file not properly formatted: invalid 'shiftreference=' line");
			}
			
			line = lineReader.readLine(); 
			matcher = libraryLine.matcher(line);
			if (matcher.matches()) {
				this.libraryFile = matcher.group(1);				
			} else {
				throw new IOException("The configuration file not properly formatted: invalid 'library=' line");
			}
			
		} finally {
			
			if (lineReader != null) {
				lineReader.close();
				
			}
		}
		
	}
	
	public void printAtoms(){
		int i,j;
		int epSize,assignmentIndex;
		System.out.println("===All atoms of the experiment including all spectra===");

		System.out.println("Index: Atoms of the whole experiment (spectrumName.expected peak) --- *:assigned");
		for (i = 0; i < this.atoms.size(); i++){
			System.out.print(this.atoms.get(i).getIndex()+": "+this.atoms.get(i).getName());
			if (this.atoms.get(i).getAssigned(this.assignment.getAssignmentWhole())== true)
				System.out.print("*");
			System.out.print(" (");
			epSize = this.atoms.get(i).getExpectedPeakSize();
			for (j = 0; j < epSize; j++){
				System.out.print(this.atoms.get(i).getExpectedPeak(j).getSpectrumName()+"."+this.atoms.get(i).getExpectedPeak(j).getIndex());
				assignmentIndex = this.atoms.get(i).getExpectedPeak(j).getAssignmentIndex();
				if (this.assignment.getAssignment(assignmentIndex)> -1)
					System.out.print("*");
				System.out.print(" ");
			}
			
			System.out.println(")");
		}
		System.out.println();
	}
	

	
	private void printConfigurationFile(){
		
		int i;
		System.out.println("*****configuration file******");
		System.out.print("peaks:=");
		for (i = 0; i < spectrumNames.length; i++){
			System.out.print(spectrumNames[i]);
			if (i!=spectrumNames.length -1)
				System.out.print(",");
		}
		System.out.println();
		System.out.println("sequence:="+this.sequenceFile);
		System.out.println("shiftreference:="+this.shiftReferenceFile);
		System.out.println("library:="+this.libraryFile);
		System.out.println("*******************");
	}

	@Override
	public String getName() {
		return "Chemical Shift Assignment Problem";
	}

	@Override
	public int getNumberOfConstraints() {
		return nConstraints;
	}

	@Override
	public int getNumberOfObjectives() {
		return nObjectives;
	}

	@Override
	public int getNumberOfVariables() {
		return nVariables;
	}


	public int getNrSpecTypes(){
		return this.numberOfSpectrumTypes;
	}
	
	public void setNrSpecTypes(int value){
		this.numberOfSpectrumTypes = value;
	}
	
	@Override
	public void close() {
		//do nothing
	}


}
