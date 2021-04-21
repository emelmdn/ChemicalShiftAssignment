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

import java.io.BufferedWriter;
import java.io.File;
import java.util.*; 
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.io.Reader;
import java.io.UnsupportedEncodingException;
import java.io.Writer;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.moeaframework.core.variable.BinaryVariable;
import org.moeaframework.util.Vector;
import org.moeaframework.util.io.CommentedLineReader;

/*
 * Public class Spectrum that saves all spectrum related attributes and their functions for all spectra of NMR.
 */
public class Spectrum{

	private ArrayList<Peak> measuredPeaks;
	private ArrayList<Peak> expectedPeaks; //Expected Peak type to be defined
	private int[] assignment;

	private Integer degeneracyValues[];
	private ArrayList<Atom> atoms;
	private int nObjectives;
	private int nConstraints;
	private float shiftNormality;
	private float alignment;
	private float degeneracy;
	private int index; //shows the index of the spectrum in ChemicalShiftAssignments.allSpectra. First one is 0.
	
	private int dimension;
	private String format;
	private ArrayList<String> identifiers; //atom names from measured peak list (#INAME lines)
	private ArrayList<String> measuredPeakIdentifiers; //atom names from measured peak list (#INAME lines)
	private ArrayList<String> expectedPeakIdentifiers; //atom names from expected peak list (#INAME lines)
	private ArrayList<Integer> indexMapper; //atom's index in a measured peak chemical shift list (0. item is 1 means, the first chemical shift of expected peak's atom (identifier) is the second one in the measured peaks 
	private String name;
	private Library library;
 	private ArrayList<String> aminoacidSequence;
 	private double[][] chemicalShiftValues; //chemical shift values of the measured peaks. [#dim][#measured peaks] atom order in the dim is the same as the one in the expected peaks
 	private int[][] sortedPeaks; //sorted measured peaks. [#dim][#measuredPeaks]. chemicalShiftValues values are sorted in quicksort2 function and this matrix show the index of chemicalShiftValues values from min to max for each dim.
 	private Tree tree;
	private int[] assignedPeakCount;  //number of expected peaks assigned to measured peaks (size = #measured peaks)
	private double d = -1; 
	private double expectedPeakCount; //expected number of measured peaks in the spectrum
	private double degeneracyFactor;
	private double pValue;
	private double[][] minumumMaximumFrequency; //min and maximum values for frequency
 	private int[] numberOfAtoms; //number of different atoms in expected peaks for dim
 	private int category = 0; //category for assignment
 	private int assignmentIndexFirstExpPeak; //index of the first expected peak of this spectrum in ChemicalShiftAssignment.assignment (First one is zero) = assignmentIndex of the first expPeak
 	private int assignmentIndexLastExpPeak; //index of the last expected peak of this spectrum in ChemicalShiftAssignment.assignment (First one is zero) = assignmentIndex of the last expPeak
	
	public Spectrum() {
		this.identifiers = new ArrayList<String>();
		this.measuredPeakIdentifiers = new ArrayList<String>();
		this.expectedPeakIdentifiers = new ArrayList<String>();
		this.indexMapper = new ArrayList<Integer>();
		this.measuredPeaks = new ArrayList<Peak>();
		this.expectedPeaks = new ArrayList<Peak>();
		this.atoms = new ArrayList<Atom>();
		this.aminoacidSequence = new ArrayList<String>();
		this.tree = new Tree();
		this.measuredPeakIdentifiers = new ArrayList<String>();
	}
	
	public Spectrum(String spectrum, int nObj, int nConst ) throws IOException {
		
		this.identifiers = new ArrayList<String>();
		this.measuredPeakIdentifiers = new ArrayList<String>();
		this.expectedPeakIdentifiers = new ArrayList<String>();
		this.indexMapper = new ArrayList<Integer>();
		this.measuredPeaks = new ArrayList<Peak>();
		this.expectedPeaks = new ArrayList<Peak>();
		this.atoms = new ArrayList<Atom>();
		this.aminoacidSequence = new ArrayList<String>();
		this.tree = new Tree();
		
		loadPeaks (spectrum+".peaks",0); //HNCA.peaks
		loadPeaks (spectrum+"_exp.peaks",1); //HNCA_exp.peaks
		
		assignPeakIndices();
		createAtoms();
		printAtoms();
		setInitialValues();
	
		
		this.nObjectives = nObj;
		this.nConstraints = nConst;
		
	}
	
	public Spectrum(String spectrum, int index) throws IOException {
		this.identifiers = new ArrayList<String>();
		this.measuredPeakIdentifiers = new ArrayList<String>();
		this.expectedPeakIdentifiers = new ArrayList<String>();
		this.indexMapper = new ArrayList<Integer>();
		this.measuredPeaks = new ArrayList<Peak>();
		this.expectedPeaks = new ArrayList<Peak>();
		this.atoms = new ArrayList<Atom>();
		this.aminoacidSequence = new ArrayList<String>();
		this.index = index;

		this.tree = new Tree();
		
		loadPeaks (spectrum+".peaks",0); //HNCA.peaks
		loadPeaks (spectrum+"_exp.peaks",1); //HNCA_exp.peaks
		assignPeakIndices();
		createAtoms();
		setInitialValues();
		
	}

	/*
	 * This function initializes  the spectrum specific values such as accuracy
	 */
	public void setInitialValues() {
		int i;
		
		this.allocateMemory();
		this.initializeValues();
		setAccuracyValues();
		initializeMaxScore();
		this.assignedPeakCount = new int[this.getMeasuredPeaksSize()];
		for (i = 0; i < this.assignedPeakCount.length; i ++){
			this.assignedPeakCount[i] = 0;
		}
		
		constructTree();
		this.sortPeaks();
		
	}
	
	/*
	 * This function initializes the sortedPeaks values for the corresponding spectrum. 
	 */
	private void initializeValues(){
		int idim, i;
		
		for (idim = 0; idim < this.dimension; idim ++){
			for (i = 0; i < this.getMeasuredPeaksSize(); i++){
				this.sortedPeaks[idim][i] = i;
			}
		}
	}
	/*
	 * This function allocates memory fo rhte peaks and the atoms of the spectrum.
	 */
	private void allocateMemory(){
		this.chemicalShiftValues = new double[this.dimension][this.getMeasuredPeaksSize()];
		this.sortedPeaks = new int[this.dimension][this.getMeasuredPeaksSize()];	
		this.numberOfAtoms = new int[this.dimension];
	}

	public double getDValue(){
		return this.d;
	}
	public void setAssignedPeakCount(int index, int value){
		this.assignedPeakCount[index] = value;
	}
	
	public void setNumberOfAtoms(int index, int value){
		if (index >= this.dimension)
			System.out.println("ERROR in setNrAtoms in Spectrum.java: Given index is bigger than the spectrum dimension!");
		this.numberOfAtoms[index] = value;
	}
	
	public int getNumberOfAtoms(int index){
		return this.numberOfAtoms[index];
	}
	
	/*
	 * This function builds the search tree to traverse the search space.
	 */
	private void constructTree(){
		int i,j, measId;
		int[] pos;
		double [][] chemicalShiftValuesTemp;

	
		//the values in the chemicalShiftValues (dimension) fits with the atom names in the expected peak atoms (0.atom of expected peak = 0. atom of the array )
		for (i = 0; i < this.getMeasuredPeaksSize(); i++){
			for (j = 0; j < this.dimension; j++){
				measId = this.indexMapper.get(j);
				this.chemicalShiftValues[j][i] = this.measuredPeaks.get(i).getChemicalShift(measId);
			}
		}
		
	
		pos = new int[this.getMeasuredPeaksSize()];
		chemicalShiftValuesTemp = new double[this.dimension][this.getMeasuredPeaksSize()];
		for (i = 0; i < this.dimension; i++){
			for (j = 0; j < this.getMeasuredPeaksSize(); j++){
				chemicalShiftValuesTemp[i][j] = this.chemicalShiftValues[i][j];
			}
		}
		
		for (i = 0; i < pos.length; i++)
			pos[i] = i;
		
		initializeBranch(this.tree, chemicalShiftValuesTemp, pos, this.dimension);
	}
	
	/*
	 * This function initializes the branches and the leaves of the  search tree.
	 */
	private void initializeBranch(Tree tree, double[][] chemicalShifts, int[] positions, int d) {

		int dimension; 
		
		dimension = this.dimension - d;
		
		chemicalShifts = quickSort(chemicalShifts, positions, dimension, 0, positions.length-1);
		
		tree.setMeasuredPeakNumber((positions.length));
		tree.setDimensions(chemicalShifts[dimension]);
		
		
		if (d == 1){
			tree.setOrientationTo(positions);
		} else if (positions.length > 1){
			tree = createAssociation(tree, chemicalShifts, positions, d, 0, (positions.length-1)); 
			
		} else{
			
			initializeLeaf(tree.getFurtherDimension(0), getColumnFromMatrix(chemicalShifts, 0), positions[0], d-1);
		
		}
		
	}
	
	public void clearEqualPeaks() {
		int i;
		for (i = 0; i < this.getExpectedPeaksSize(); i++) {
			this.getExpectedPeak(i).clearEqualPeaks();
		}
	}
	

	/*
	 * This function creates the relevant associates in the tree among its nodes. 
	 */
	private Tree createAssociation(Tree tree, double[][] chemShifts, int[] positions, int d, int begin, int end) {
		int medium;
		double [][] newChemShift;
		int row;
		int [] newPosition;
		
		medium = (begin + end)/2;

		row = chemShifts.length;
		newChemShift = new double[row][end-begin+1];
		newPosition = new int[end-begin +1];
	
		newChemShift = copyColumnsToMatrixFromMatrix(chemShifts, begin, end);
		newPosition = copyValuesToArrayFromArray(positions, begin, end);
		
		
		initializeBranch(tree.getFollowingDimension(medium), newChemShift, newPosition, d-1);
					
		if (begin < medium){
			tree = createAssociation(tree, chemShifts, positions, d, begin, medium);
		} else{
			initializeLeaf (tree.getFurtherDimension(medium),getColumnFromMatrix(chemShifts, medium), positions[medium], d-1);
		}
		
		if ((medium+1) < end){
			tree = createAssociation(tree, chemShifts, positions, d, medium+1, end);
		} else{
			initializeLeaf(tree.getFurtherDimension(medium+1), getColumnFromMatrix(chemShifts, (medium+1)), positions[medium+1], (d-1));
		}
			
		return tree;
	}

	public Tree getTree(){
		
		return this.tree;
	}
	
	/*
	 * This function adds one leaf to the existing tree.
	 */
	private void initializeLeaf(Tree r, double[] ppm, int pos, int d) {
		int dimension;
		dimension = ppm.length -d;
		
		r.setFollowingDimensionToSingel(ppm[dimension]);
		
	
		if (d == 1){
			r.addToOrientation(pos);
		} else {
			initializeLeaf(r.getFurtherDimension(0), ppm, pos, d-1);
		}
		
	
	}

	/*
	 * The values of the position array from[start] to [end] are copied into a new array and sent as result.
	 */
	private int[] copyValuesToArrayFromArray(int[] position, int start, int end){
		int[] newPosition;
		int newSize = end - start + 1;
		int i;
		
		newPosition = new int[newSize];
		for (i = 0; i < newSize; i++)
			newPosition[i]= position[start+i];
		return newPosition;
	}
	/*
	 * This funciton copies given column of one matrix to another and sends the new matrix as result.
	 * newppm = chemShift(:, start,end)
	 */
	private double[][] copyColumnsToMatrixFromMatrix(double[][] chemShift, int start, int end) {
		int i,j, k;
		int row;
		int newSize = end - start + 1;
		double[][] newChemShift;
		
		
		row = chemShift.length;
	
		newChemShift = new double[row][newSize];
		
		k = 0;
		for (i = start; i <= end; i++){
			for (j = 0; j < row; j++){	
				newChemShift[j][k] = chemShift[j][i];
				
			}
			k++;
		}
		
		
		return newChemShift;
	}

	private double[][] quickSort(double[][] chemShiftArray, int[] positions, int dimension, int start, int end){
		
		int i,j, h2;
		double [] chemShifts;
		double med;
		double[] h1;
	
		i = start;
		j = end;
		
		chemShifts = new double[chemShiftArray[0].length];
		chemShifts = chemShiftArray[dimension];
		med = chemShifts[(start + end)/2];
		h1 = new double[chemShiftArray.length];
		while (true){
			while (chemShifts[i] < med){
				i++;
			}
			while (med < chemShifts[j]){
				j--;
			}
			if (i <= j){
				h1 = getColumnFromMatrix(chemShiftArray,i);
				chemShiftArray = copyColumnToMatrixFromMatrix(chemShiftArray,i,chemShiftArray,j);
				chemShiftArray = copyColumnToMatrixFromArray(chemShiftArray, j, h1);
				h2 = positions[i];
				positions[i] = positions[j];
				positions[j] = h2;
				i++;
				j--;
			}
			if (i > j){
				break;
			}
		}
		
		if (start < j){
			chemShiftArray = quickSort(chemShiftArray, positions, dimension, start, j);
		}
		
		if (i < end){
			chemShiftArray = quickSort(chemShiftArray, positions, dimension, i, end);
		} 
		return chemShiftArray;
	}

	/*
	 * This function copies FROM j. column of source TO i. column of destination 
	 * 
	 * 
	 */
	private double[][] copyColumnToMatrixFromArray(double[][] destination, int j, double[] source) {
		int a;
		int row = destination.length;
		
		for (a = 0; a < row; a++){
			destination[a][j]= source[a];
		}
		return destination;
		
	}

	/*
	 * 
	 * This function copies FROM j. column of source TO i. column of destination 
	 */
	private double[][] copyColumnToMatrixFromMatrix(double[][] destination, int i, double[][] source, int j) {
		int a, b;
		int row = source.length;
		
		for (a = 0; a < row; a++){
			destination[a][i]= source[a][j];
		}
		return destination;
		
	}
	/*
	 * This function returns i. column of the matrix
	 * 
	 */
	private double[] getColumnFromMatrix(double[][] matrix, int i) {
		int a;
		int row = matrix.length;
	
		double[] result = new double[row];
		for (a = 0; a < row; a++)
			result[a] = matrix[a][i];
		
		return result;
	}

	/*
	 * This function calculates maxScore, which is denominator of the score formula.
	 * It is assignment independent and is a constant value during whole calculation (as soon as expected peaks do not change)
	 * This function is used only for a single spectrum experiment. In case of multiple spectrum, new one in chemicalShiftAssignment class is applied
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
	}

	/*
	 * This function calculates the degeneracy values for the given assignment.
	 * Degeneracy is set to -1 if not assigned, o/w shows the total number of assignment
	 * No need to initialize the degeneracy values before because it resets the deg values first.
	 */
	public void calculateDegeneracy(int[] asg){
		int degSize = asg.length;
		boolean calculated[] = new boolean[degSize];
		int i, j, count;
		
		this.degeneracyValues = new Integer[degSize];
		
		for (i = 0; i < degSize; i++){ // assign initial values
			calculated[i] = false;
			this.degeneracyValues[i] = -1;
		}
		
		for (i = 0; i < degSize; i++){
			if (calculated[i] == false){
				if (asg[i] == -1){ // the expected peak is  not assigned
					calculated[i] = true;
				} else {
					count = 0;
					for (j = 0; j < degSize; j++){
						if (asg[i] == asg[j]){
							count ++;
						}
					}
					for (j = 0; j < degSize; j++) {
						if (asg[j] == asg[i]){
							calculated[j] = true;
							this.degeneracyValues[j] = count;
						}
					}
				}
			}
		}
		
	}
	
	/*
	 * 
	 * This funciton sets the number of objectives and constraints as given. 
	 */
	public void setObjectiveConstraintCounts(int objCount, int conCount){
		this.nObjectives = objCount;
		this.nConstraints = conCount;
	}
	
	/*
	 * This function goes through all assigned values of the atom and finds the real average of these chemical frequencies
	 * IF change: apply to calculateAverageFrequency, calculateG
	 */
	public void calculateAverageFrequencyOfAtom(Atom at, int[] asg, int i){
		
		int j, epSize, expPeakIndex, measuredPeakIndex, dimIndex; 
		epSize = at.getExpectedPeakSize();
		Peak expPeak = new Peak();
		Peak measuredPeak = new Peak();
		double currentFrequency;
		
		for (j = 0; j < epSize; j++){ //for all expected peaks of the atom
			expPeak = at.getExpectedPeak(j);
			if (expPeak.getAssigned() == true){
				expPeakIndex = expPeak.getIndex(); 
				measuredPeakIndex = asg[expPeakIndex];	
				measuredPeak = this.measuredPeaks.get(measuredPeakIndex);
				dimIndex = at.getExpectedPeaksIndexMapperOf(j);
				currentFrequency = measuredPeak.getChemicalShift(dimIndex);
				at.addMeasuredFrequency(currentFrequency);
			}
		}
		at.calculateAverageFrequency();
	
	}
	
	/*
	 * This function is used during the Global score calculation only for one spectrum case.
	 */
	public double cdf(double x, double d){
		double cdf, e;
		Derf a = new Derf();
				
		e = x/(Math.sqrt(2.0)*d);
		cdf = 1.0 - 0.5*(a.erfcc(e));
		
		return cdf;
	}
	

	/*
	 * This function resets all of the atom assignments. 
	 */
	public void resetAtomAssignments(){
		int i,j;
		int atomCount;
		int epCount;
		
		atomCount = this.atoms.size(); 
		
		for (i = 0; i < atomCount; i++){
			this.atoms.get(i).setAssigned(false);
			epCount = this.atoms.get(i).getExpectedPeakSize();
			for (j = 0; j < epCount; j++){
				this.atoms.get(i).getExpectedPeak(j).setAssigned(false);
			}
		}
	}

	
	/*
	 * This function assigns punishment for the constraints and returns them as result.
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
	

	
	/* This function loads measured peaks and expected peaks from the fileName
	 * It reads either expected or measured peaks depending on type input parameter: 0 ==> measured peak, 1 ==> expected peak
	 */
	public void loadPeaks(String fileName, int type) throws IOException{
		
		int i, dimension, peakNumber, colorIndex, colorCode, anInt; 
		double[] chemicalShifts;
		String[] atoms;
 		double peakVolume, peakVolumeError;
		String spectrumType,integrationMethod;
		Peak mp;
		Pattern dimensionLine = Pattern.compile("# Number of dimensions (\\d+)");
		Pattern formatLine = Pattern.compile("#FORMAT (.*?)");
		Pattern identifierLine = Pattern.compile("#INAME (\\d+)\\s(.*?)");
		Pattern spectrumLine = Pattern.compile("#SPECTRUM (.*?) (.*?)");
		
		CommentedLineReader lineReader = null;
		String line = null;
		Matcher matcher = null;
		
		try{
			File file = new File(fileName);
			if (file.exists() == false){
				throw new IOException("Input peak file "+fileName+" does not exist!");
			}
			Reader reader = new FileReader(file);
			lineReader = new CommentedLineReader(reader);
			line = lineReader.readLine(); //the number of dimension line
			matcher = dimensionLine.matcher(line);
			if (matcher.matches()) {
				setDimension(Integer.parseInt(matcher.group(1)));		
				dimension = getDimension();
			} else {
				throw new IOException("XEASY peak list file "+fileName+" not properly formatted: invalid '# Number of dimension' line");
			}
			chemicalShifts = new double[dimension];
			atoms = new String[dimension];
		
			line = lineReader.readLine(); 
			matcher = formatLine.matcher(line);
			if (matcher.matches()) {
				setFormat(matcher.group(1));	
			} else {
				throw new IOException("XEASY peak list file not properly formatted: invalid 'format' line");
			}
			
			//identifier lines
			for (i = 0; i < dimension; i++){
				line = lineReader.readLine(); 
				matcher = identifierLine.matcher(line);
				if (matcher.matches()) {
					if (type == 0)//measured peak
						addToMeasuredPeakIdentifiers(matcher.group(2));
					if (type == 1) //expected peak
						addToExpectedPeakIdentifiers(matcher.group(2));
					addToIdentifiers(matcher.group(2));
				} else {
					throw new IOException("XEASY peak list file not properly formatted: invalid 'identifier' line");
				}
			}
			
			//spectrum line
			line = lineReader.readLine(); 
			matcher = spectrumLine.matcher(line);
			if (matcher.matches()) {
				setName(matcher.group(1));	
			} else {
				throw new IOException("XEASY peak list file not properly formatted: invalid 'spectrum' line");
			}
			
			line = lineReader.readLine();
			String peakLinePattern= "(\\s)*(\\d+)(\\s)*";
			for (i=0; i<dimension; i++){
				peakLinePattern += "(-?\\d+\\.\\d+)(\\s)*";
			}
			peakLinePattern += "(\\d+)(\\s)*(\\w)(\\s)*((-?\\d+(\\.\\d+)?)?(E))([-+]?\\d+)?(\\s)*((-?\\d+(\\.\\d+)?)?(E))([-+]?\\d+)?(\\s)*";
			peakLinePattern += "(\\s)*(\\w)(\\s)*(\\d+)"; 
			for (i =0; i<dimension; i++){
				peakLinePattern += "( |\t)+(\\S+)*";
			}
			peakLinePattern += "(.*?)";
			
			Pattern peakLine = Pattern.compile(peakLinePattern);
			
			while(line!=null){
				matcher = peakLine.matcher(line);
				if (matcher.matches()) {
					peakNumber = Integer.parseInt(matcher.group(2));
					for (i = 0 ; i < dimension; i++){
						chemicalShifts[i] = Double.parseDouble((matcher.group(4 + (i*2))));
					}
				} else {
					throw new IOException("XEASY peak list file ["+fileName+"] not properly formatted: invalid 'peak' line OR empty line");
				}
				colorIndex = 4 + (dimension*2);
				colorCode = Integer.parseInt(matcher.group(colorIndex));
				
				spectrumType = matcher.group(colorIndex + 2);
				peakVolume = Double.parseDouble(matcher.group(colorIndex + 5)) * (Math.pow(10, Integer.parseInt(matcher.group(colorIndex + 8))));
				peakVolumeError = Double.parseDouble(matcher.group(colorIndex + 11)) * (Math.pow(10, Integer.parseInt(matcher.group(colorIndex + 14))));

				integrationMethod = matcher.group(colorIndex + 17);
				
				anInt = Integer.parseInt(matcher.group(colorIndex+19));
				
				String s;
				for (i = 0 ; i < dimension; i++){
					s = matcher.group(colorIndex + 21 + (i*2));
					atoms[i] = s;
					
				}
				mp = new Peak(dimension, peakNumber,chemicalShifts, colorCode, spectrumType, peakVolume, peakVolumeError,integrationMethod, anInt, atoms, this.getName(), this.getIndex());
				if (type == 0){ //measured peak
					this.addMeasuredPeak(mp);
				} else if (type == 1){ //expected peak
					this.addExpectedPeak(mp);
				} else{
					throw new IOException("Undefined input type while reading the spectrum input files");
				}
				mp = null;
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
	 * This function assigns indices of measured peaks
	 */
	public void assignPeakIndices(){
		int i, mpSize, epSize;
		
		mpSize = this.getMeasuredPeaksSize();
		epSize = this.getExpectedPeaksSize();
		
		//asssign indices of measured peaks
		for (i = 0; i < mpSize; i++){
			this.measuredPeaks.get(i).setIndex(i);
		}
		
		for (i = 0; i < epSize; i++){
			this.expectedPeaks.get(i).setIndex(i);
		}
	}
	
	/*
	 * This function creates the atoms of the spectrum and assigns the expected peaks to them.
	 * This is static for the spectrum. Assignment affects are not evaluated yet.
	 */
	public void createAtoms() throws IOException{
		
		int i, j, atomIndex;
		String atomName, genericAtomName; 
		Peak pk;
		Atom at; 
		int residueNumber;

		try{
			createIndexMapper();		
			for (i = 0; i < this.expectedPeaks.size(); i++){
				pk = new Peak(this.expectedPeaks.get(i));
				for (j = 0; j < this.dimension; j++){
					atomName = this.expectedPeaks.get(i).getAtom(j).getName();
					residueNumber = Integer.parseInt(atomName.substring(atomName.indexOf(".")+1));
					genericAtomName = atomName.substring(0, atomName.indexOf("."));
					
					atomIndex = atomExists(atomName);
					if (atomIndex == -1){ //atom is new
						this.expectedPeaks.get(i).getAtom(j).setIndex(this.atoms.size());
						at = new Atom(atomName, indexMapper.get(j), residueNumber, genericAtomName);
						at.setIndex(this.atoms.size());
						at.addExpectedPeak(pk,indexMapper.get(j));
						this.atoms.add(at);
					} else if (atomIndex >= 0){ // atom exists
						this.expectedPeaks.get(i).getAtom(j).setIndex(atomIndex);
						this.atoms.get(atomIndex).addExpectedPeak(pk,indexMapper.get(j));
					
					} else{
						throw new IOException("Atom index is not found for expected peak in createAtoms method of Spectrum.class");
					}		
				}	
			}
		} finally {
	
		}
	
	}
	
	/*
	 * This function creates the values for the index mapper from chemical shift indices of atoms in expected peak into the chemical shift on measured peak
	 * For ex: if 0. element's value is 2: 0.atom of the expected peak is at the 2. element of measured peak. 
	 * (used to get chemical shift values from measured peak for the atom)
	 */
	private void createIndexMapper() {
		int i, j;
		for (i = 0; i<this.dimension; i++){
			for (j = 0; j< this.dimension; j++){
				if (this.expectedPeakIdentifiers.get(i).equals(this.measuredPeakIdentifiers.get(j))){
					this.indexMapper.add(j);
					break;
				}
			}
		}
	}


	/*
	 * This function sets the accuracy values for each atom of the spectrum.
	 */
	public void setAccuracyValues(){
		int i, size;
		size = this.atoms.size();
		
		for (i = 0; i < size; i++){
			this.atoms.get(i).setAccuracyValues();
		}
	}
	
	/*
	 * This function returns -1 if atom does not exist in spectrum
	 * returns the index of the atom otherwise
	 */
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
	 * This function prints all atoms of the spectrum, where
	 * the assigned ones are marked with * 
	 */
	public void printAtoms(){
		int i,j;
		int epSize;
		System.out.println("Index: Atoms of the spectrum "+ this.getName()+" (their expected peaks) --- *:assigned");
		for (i = 0; i < this.atoms.size(); i++){
			System.out.print(this.atoms.get(i).getIndex()+": "+this.atoms.get(i).getName());
			if (this.atoms.get(i).getAssignedOld()== true)
				System.out.print("*");
			System.out.print(" (");
			epSize = this.atoms.get(i).getExpectedPeakSize();
			for (j = 0; j < epSize; j++){
				System.out.print(this.atoms.get(i).getExpectedPeak(j).getIndex());
				if (this.atoms.get(i).getExpectedPeak(j).getAssigned() == true)
					System.out.print("*");
				System.out.print(" ");
			}
			
			System.out.println(")");
		}
		System.out.println();
	}
	
	/*
	 * This function prints all values of the spectrum.
	 */
	public void printSpectrum(){
		
		int dimension = getDimension();
		int i, mp;
		System.out.println("************SPECTRUM:"+name+"(start)**********");
		System.out.println("# Number of dimensions "+ dimension);
		System.out.println("#FORMAT "+ format);
		for (i = 0; i < dimension; i++){
			System.out.println("#INAME "+(i+1)+" "+ identifiers.get(i));
		}
		System.out.println("#SPECTRUM " + name);
		mp = measuredPeaks.size();
		for (i= 0; i < mp; i++){
			this.measuredPeaks.get(i).printPeakLine();
		}
		System.out.println("Expected Peaks:");
		mp = expectedPeaks.size();
		for (i= 0; i < mp; i++){
			this.expectedPeaks.get(i).printPeakLine();
		}
		System.out.println("************SPECTRUM:"+name+"(end)**********");
	}
	

	public ArrayList<Peak> getMeasuredPeaks() {
		return measuredPeaks;
	}
	
	public int getMeasuredPeaksSize() {
		return measuredPeaks.size();
	}

	
	public int getAssignedPeakCount(int i){
		return this.assignedPeakCount[i];
	}
	
	
	public int[] getAllAssignedPeaks(){
		return this.assignedPeakCount;
	}
	public void increaseAssignedPeakCount(int measId){
		this.assignedPeakCount[measId] += 1;
	}

	public void setMeasuredPeaks(ArrayList<Peak> measuredPeaks) {
		this.measuredPeaks = measuredPeaks;
	}
	public ArrayList<Peak> getExpectedPeaks() {
		return expectedPeaks;
	}
	
	public Peak getExpectedPeak(int index) {
		return expectedPeaks.get(index);
	}
	
	public int getExpectedPeaksSize() {
		return expectedPeaks.size();
	}
	
	/*
	 * This function returns the size of the assignment.
	 * The assignment size must be the same as the expected peak size
	 */
	public int getAssignmentSize(){
		return expectedPeaks.size();
	}

	public double getChemicalShiftValues(int i, int j){
		return this.chemicalShiftValues[i][j];
	}
	
	public void setExpectedPeaks(ArrayList<Peak> expectedPeaks) {
		this.expectedPeaks = expectedPeaks;
	}

	public int getnObjectives() {
		return nObjectives;
	}

	public void setnObjectives(int nOfObjectives) {
		this.nObjectives = nOfObjectives;
	}

	public int getnConstraints() {
		return nConstraints;
	}

	public void setnConstraints(int nOfConstraints) {
		this.nConstraints = nOfConstraints;
	}

	public float getShiftNormality() {
		return shiftNormality;
	}

	public void setShiftNormality(float shiftNormality) {
		this.shiftNormality = shiftNormality;
	}

	public float getAlignment() {
		return alignment;
	}

	public void setAlignment(float alignment) {
		this.alignment = alignment;
	}
	
	public 	ArrayList<Atom> getAtoms(){
		return this.atoms;
	}

	public Atom getAtom(int i){
		return this.atoms.get(i);
	}
	
	public float getDegeneracy() {
		return degeneracy;
	}

	public void setDegeneracy(float degeneracy) {
		this.degeneracy = degeneracy;
	}

	public int getDimension() {
		return dimension;
	}

	public void setDimension(int dimension) {
		this.dimension = dimension;
	}

	public String getFormat() {
		return format;
	}

	public void setFormat(String format) {
		this.format = format;
	}

	public ArrayList<String> getIdentifiers() {
		return identifiers;
	}

	public String getIdentifier(int i) {
		return identifiers.get(i);
	}
	
	public void setIdentifiers(ArrayList<String> identifiers) {
		this.identifiers = identifiers;
	}
	
	public void addToExpectedPeakIdentifiers (String identifier) {
		this.expectedPeakIdentifiers.add(identifier);
	}

	public void addToMeasuredPeakIdentifiers (String identifier) {
		this.measuredPeakIdentifiers.add(identifier);
	}
	public void addToIdentifiers (String identifier) {
		this.identifiers.add(identifier);
	}
	public String getName() {
		return name;
	}

	public void setName(String name) {
		this.name = name;
	}
	
	public void addMeasuredPeak(Peak m){
		Peak mp = new Peak(m);
		this.measuredPeaks.add(mp);
	}
	
	public void addExpectedPeak(Peak m){
		Peak mp = new Peak(m);
		this.expectedPeaks.add(mp);
	}

	public Library getLibrary() {
		return library;
	}

	public void setLibrary(Library library) {
		this.library = new Library(library);
	
	}

	/*
	 * This function assigns the reference frequency values from CSTable
	 */
	public void assignReferenceFrequency(){
	
		int i, residueNumber;
		Pattern pt = Pattern.compile("(\\S+)*(\\.)(\\d+)");
		Matcher matcher;
		String atomName, residueName;
		Double ref, dev;
		
		for (i = 0; i < this.atoms.size(); i++){ // for all atoms 
			matcher = pt.matcher(this.atoms.get(i).getName());
			if (matcher.matches()){
				atomName = matcher.group(1);
				residueNumber = Integer.parseInt(matcher.group(3));
				residueName = getResidueName (residueNumber);
				ref = this.library.getReferenceFrequency(residueName, atomName);
				dev = this.library.getReferenceDev(residueName, atomName);
				this.atoms.get(i).setReferenceFrequency(ref);
				this.atoms.get(i).setCSTableFreq(dev);
			}
		}
	}

	/*
	 * This function returns the existing minimum and maximum chemical shift frequency values. 
	 */
	public void getMinMaxFreq(){
		
		int idim;
		this.minumumMaximumFrequency = new double [this.dimension][2];
		

		for (idim = 0; idim < this.dimension; idim++){
			this.minumumMaximumFrequency[idim][0] = this.chemicalShiftValues[idim][this.sortedPeaks[idim][0]];
			this.minumumMaximumFrequency[idim][1] = this.chemicalShiftValues[idim][this.sortedPeaks[idim][this.getMeasuredPeaksSize()-1]];
		}
		
	}
	
	/*
	 * This function sorts the chemical shift frequency values of the spectrum.
	 */
	private void sortPeaks(){
		
		int idim;
		int[] pos = new int[this.getMeasuredPeaksSize()];
		double[] ppm = new double[this.getMeasuredPeaksSize()];
		
		
		for (idim = 0; idim < this.dimension; idim ++){

			pos = this.sortedPeaks[idim];
			ppm =  this.chemicalShiftValues[idim];
			
			
			quickSort2(ppm, pos, 0, this.getMeasuredPeaksSize()-1);
			
		}
		
	}
	

	
	/*
	 * Quick sort algorithm. 
	 * The values are saved in the vals. They are not changed at all. Result is saved in the pos array, where the indices
	 * of the sorted values are saved.
	 * 
	 */
	public void quickSort2(double[] vals, int[] pos, int beg, int end){

		int i, j, help;
		double med;
		
		i = beg;
		j = end;
		med = vals[pos[(i+j)/2]];
		
		while (i <= j){
			while (vals[pos[i]] < med)
				i ++;
			while (med < vals[pos[j]])
				j--;
			if (i <= j){
				help = pos[i];
				pos[i] = pos[j];
				pos[j] = help;
				i++;
				j--;
			}
		}
		if (beg < j)
			quickSort2(vals, pos, beg, j);
		if (i < end)
			quickSort2(vals, pos, i, end);		
	}
	
	
	public void setAminoAcidSequence(ArrayList<String> seq) {
		int i;
		for (i = 0; i < seq.size(); i++){
			this.aminoacidSequence.add(seq.get(i));
		}
		
	}
	
	/*
	 *  This function sends the number of the residue. (1 is the first one)
	 */
	public String getResidueName(int residueNumber){
		
		
		residueNumber = residueNumber - 1;
		return this.aminoacidSequence.get(residueNumber);
	}
	
	public int getAtomSize(){
		return this.atoms.size();
	}
	
	public void printIdentifiers(){
		int i;
		
		System.out.println("size: " + this.dimension);
		System.out.println("expected identifiers size: "+ this.expectedPeakIdentifiers.size());
		System.out.println("----Expected Peak Identifiers----");
		for (i=0; i<this.dimension; i++)
			System.out.print(expectedPeakIdentifiers.get(i) + " ");
		System.out.println();

		System.out.println("----Measured Peak Identifiers----");
		for (i=0; i<this.dimension; i++)
			System.out.print(measuredPeakIdentifiers.get(i)+" ");
		System.out.println();
		

		System.out.println("----  Identifiers----");
		for (i=0; i<this.dimension; i++)
			System.out.print(identifiers.get(i)+ " ");
		System.out.println();
	}

	public int[] getAssignment() {
		return assignment;
	}

	public void setAssignment(int[] asg) {
		this.assignment = new int [asg.length];
		this.assignment = asg;
	}

	public ArrayList<Integer>  getIndexMapper() {
		return this.indexMapper;
	}

	/*
	 * This function returns the index of the atom in the expected peak of the spectrum. for ex: HN atom's index in HNCA spectrum is 0 (H N CA)
	 * This function is used to find parameter for indexMapper logic for spectrum extend 
	 */
	
	public int findExpectedPeakIndexOfAtom(String atomName, Peak expectedPeak) {
		int j;
		String atomNameInExpPeak;
		
		atomName =  removeResidueIndexFromAtomName(atomName);
		System.out.println("Search atom "+ atomName);
		for (j = 0; j < this.dimension; j++){
			atomNameInExpPeak = removeResidueIndexFromAtomName(this.expectedPeaks.get(0).getAtom(j).getName());
			System.out.println(atomName+" =? "+ atomNameInExpPeak);
			if (atomName.equals(atomNameInExpPeak)){ //if the atom is found 
				System.out.println("Same!");
				return j;			//return its index
			}else
				System.out.println("Not same!");
		}

		return -1; //if not found, return error
	}
	
	/*
	 * This function removes the last and the integer afterwards which shows the residue number so that atom names can be compared
	 * for ex: HA.1 is changed to HA 
	 */
	
	public String removeResidueIndexFromAtomName(String atomName){
		
		if (atomName.contains(".")){ //if atom name contains dot.
			atomName =  atomName.substring(0, atomName.lastIndexOf(".")); //remove the values after dot.
		}
		return atomName;
	}

	/*
	 * This function is for extend spectrum. It copies all expected peaks of fromAtom to toAtom
	 * This is to copy all expected peaks of a spectrum to the chemicalshiftAssignment.atoms.expectedPeaks
	 */
	public void copyExpectedPeaksFromOneAtomToAnotherAtom(Atom fromAtom, Atom toAtom, int startIndex) {
		
		int i;
		for (i = 0; i < fromAtom.getExpectedPeakSize(); i++){
			toAtom.addExpectedPeakWithAssignmentIndex(fromAtom.getExpectedPeak(i), fromAtom.getExpectedPeaksIndexMapperOf(i), startIndex);
		}
		
	}

	public int getIndex() {
		return index;
	}

	public void setIndex(int index) {
		this.index = index;
	}
	
	public Integer[] getDegeneracyValues(){
		return this.degeneracyValues;
	}
	
	public int getDegeneracyValue(int i){
		return (int)this.degeneracyValues[i];
	}

	public void setpExistSum(double a){
		this.expectedPeakCount = a;
	}
	
	public double getExpectedPeakCount(){
		return this.expectedPeakCount;
	}
	public void setDegeneracyFactor(double a){
		this.degeneracyFactor = a;
	}
	
	public double getDegeneracyFactor(){
		return this.degeneracyFactor;
	}
	
	public void setPValue(double a){
		this.pValue = a;
	}
	
	public double getPValue(){
		return this.pValue;
	}
	
	public int getCategory(){
		return this.category;
	}
	
	public void setCategory (int value){
		this.category = value;
	}

	public int getAssignmentIndexFirstExpPeak() {
		return assignmentIndexFirstExpPeak;
	}

	public void setAssignmentIndexFirstExpPeak(int assignmentIndexFirstExpPeak) {
		this.assignmentIndexFirstExpPeak = assignmentIndexFirstExpPeak;
	}

	public int getAssignmentIndexLastExpPeak() {
		return assignmentIndexLastExpPeak;
	}

	public void setAssignmentIndexLastExpPeak(int assignmentIndexLastExpPeak) {
		this.assignmentIndexLastExpPeak = assignmentIndexLastExpPeak;
	}

}
