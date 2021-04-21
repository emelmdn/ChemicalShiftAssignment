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

import java.io.IOException;
import java.math.BigDecimal;
import java.sql.PseudoColumnUsage;
import java.util.ArrayList;

//import sun.misc.FloatingDecimal;

/*
 * Public class to save all atom related attributes and implement their functions
 */
public class Atom {
	private int fileNumber;
	private int index; //index of the atom in the Spectrum.atoms ArrayList, first atom's index is zero
	private String name;
	private String type;
	private boolean assigned; 
	private float[] coordinates = new float[3]; // x-, y- and z-coordinates 
	private int[] covalentConnectivityAtoms = new int[4]; //	four atom numbers indicating covalent connectivities (if there are less than four connectivities, the corresponding numbers are set to 0)
	private int pseudoAtomNumber; //	the atom number of the corresponding pseudo atom (or 0 if there is no corresponding pseudo atom)
	private ArrayList<Peak> expectedPeaks;
	private ArrayList<Integer> expectedPeaksIndexMapper;// the order of the atom in the measured peak(for ex: 0 means that this atom is 0. value of the chemical shifts of measured peak)
	private ArrayList<Double> measuredFrequencies;
	private int measuredPeakIndex; // the order of the atom in the measured peak(for ex: 0 means that this atom is 0. value of the chemical shifts of measured peak)
	private double averageFrequency = 0;
	private double accuracy; // will be set in the setAccuracy method depending on the atom name
	private double deviationValue = 0;
	private double frequencyValue = 0;
	private double referenceFrequency = 0; //reference frequency of the atom from CSTABLE  
	private double csTableFrequency = 0; // dev from CSTable. It is the one next to frequency
	private double globalScore = 0;
	private String atomIdentifier; //HN, N or C for HNCA spectrum. The atom identifier name of the atom in the spectrum definition
	private int atomsIndex; //The index of this atom in the ChemicalShiftAssignment.atoms
	private int residueNumber; //The residue number of the atom. The first one is 1
	private String genericAtomName; //generic atom name for ex: H, CA, ...etc
	private ArrayList<Integer> indices; //the indices of that atom in each spectrum's atoms array. For ex: if the first value is 3 (let us assume first spectrum is HNCA) then this atom's index in HNCA spectrum.atoms is 3.
	private double lowSearchBoundary;
	private double highSearchBoundary;
	private int assignedMeasuredPeakCount = 0;  //number of assigned measured peaks to the atom in the whole landscape
	private double measuredPeakFrequencySum = 0;  //sum of assigned measured peaks' frequencies to the atom in the whole landscape
	private double measuredPeakAverage = 0;  //average of assigned measured peaks' frequencies to the atom in the whole landscape
	private int assignedMeasuredPeakCountSpectrum = 0; //number of assigned measured peaks to the atom in the specific spectrum
	private double measuredPeakFrequencySumSpectrum = 0; //sum of assigned measured peaks' frequencies to the atom in the specific spectrum
	private double measuredPeakAverageSpectrum = 0;  //average of assigned measured peaks' frequencies to the atom in the specific spectrum
	private double dynamicExpectedFrequencyDeviation = 0;
	private double dynamicExpectedFrequency = 0;
	private double dynamicMeasuredFrequencyValue = 0;
	private double dynamicMeasuredFrequency = 0;
	private boolean scoreAssigned; // used in updateLocalScore only
	private double localScoreValue = 0 ; // used in updateLocalScore only 
	private int phaseNr = 0;  // used in updateLocalScore only 
	private ArrayList<Integer> expectedPeakListAssigned; // used in updateLocalScore only 
	private double averageScoreLocal = 0.5; // used in updateLocalScore only 
	private double deviationLocal = 0.2; // used in updateLocalScore only 
	private double deviationMeanValue = 1.0;//used only in updateLocalScore
	private double maximumScore = 0;//used only in updateLocalScore
	private int xyz = 0;
	private double probabilitySum = 0; //sum of prob of all exp Peaks in all spec
	private double[] probabilitySumSpectrum; // sum of prob of exp Peaks in each spectrum
	private double[] scalingFactorSpectrum; //for test of weighting of spectra, scaling factor
			
	public Atom() {
		this.assigned = false;
		this.expectedPeaks = new ArrayList<Peak>();
		this.expectedPeaksIndexMapper = new ArrayList<Integer>();
		this.measuredFrequencies = new ArrayList<Double>();
		this.indices = new ArrayList<Integer>();
	}

	/*
	 * Constructor with the given parameters.
	 */
	public Atom(int fileNumber,	String name, String type, float[] coordinates, int[] covalentConnectivityAtoms, int pseudoAtomNumber){
	
		int i;
		this.assigned = false;
		this.expectedPeaks = new ArrayList<Peak>();
		this.expectedPeakListAssigned = new ArrayList<Integer>();
		this.expectedPeaksIndexMapper = new ArrayList<Integer>();
		this.measuredFrequencies = new ArrayList<Double>();
		this.fileNumber =  fileNumber;
		this.name = name;
		this.type = type;
		for (i = 0; i < 3; i++)
			this.coordinates[i] = coordinates[i];
		for (i = 0; i < 4; i++)
			this.covalentConnectivityAtoms[i] = covalentConnectivityAtoms[i];
		this.pseudoAtomNumber = pseudoAtomNumber;		
		this.indices = new ArrayList<Integer>();
	}
	
	/*
	 * Copy constructor
	 */
	public Atom (Atom at){
		
		int i;
		
		this.assigned = false;
		this.expectedPeaks = new ArrayList<Peak>();
		this.expectedPeakListAssigned = new ArrayList<Integer>();
		this.expectedPeaksIndexMapper = new ArrayList<Integer>();
		this.measuredFrequencies = new ArrayList<Double>();
		this.indices = new ArrayList<Integer>();
		this.fileNumber = at.fileNumber;
		this.index = at.index;
		this.name = at.name;
		this.type = at.type;
		this.assigned = at.assigned;
		for (i = 0; i < 3; i++)
			this.coordinates[i] = at.coordinates[i];
		for (i = 0; i < 4; i++)
			this.covalentConnectivityAtoms = at.covalentConnectivityAtoms;
		this.pseudoAtomNumber = at.pseudoAtomNumber;
		for (i = 0; i < at.getExpectedPeakSize(); i++){
			Peak expectedPeak = new Peak(at.getExpectedPeak(i));
			this.expectedPeaks.add(expectedPeak);
		}
		for (i = 0; i < at.expectedPeaksIndexMapper.size(); i++){
			this.expectedPeaksIndexMapper.add(at.getExpectedPeaksIndexMapperOf(i));
		}
		for (i = 0; i < at.getMeasuredFrequencySize(); i++){
			this.addMeasuredFrequency(at.getMeasuredFrequency(i));
		}
		this.averageFrequency = at.averageFrequency;
		this.accuracy = at.accuracy;
		this.deviationValue = at.deviationValue;
		this.referenceFrequency = at.referenceFrequency;
		this.frequencyValue = at.frequencyValue;
		this.atomIdentifier = at.atomIdentifier;
		this.measuredPeakIndex = at.measuredPeakIndex;
		this.atomsIndex = at.atomsIndex;
		this.residueNumber = at.residueNumber;
		this.genericAtomName = at.genericAtomName;
		for (i = 0; i < at.indices.size(); i++){
			this.indices.add(at.getIndicesSingleSpectrum(i));
		}
		this.lowSearchBoundary = at.lowSearchBoundary;
		this.highSearchBoundary = at.highSearchBoundary;
	}
	
	/*
	 * Constructor with the given parameters.
	 */
	public Atom(String name, int spectrumCount, int spectrumIndex, int a, int residueNumber, String genericAtomName){
		this.name = name;
		this.indices = new ArrayList<Integer>();
		this.expectedPeakListAssigned = new ArrayList<Integer>();

		for (int i = 0; i < spectrumCount; i++)
			indices.add(-1);
		indices.set(spectrumIndex, a);	
		this.index = a;
		this.residueNumber = residueNumber;
		this.genericAtomName = genericAtomName;
		this.expectedPeaks = new ArrayList<Peak>();
		this.expectedPeaksIndexMapper = new ArrayList<Integer>();
		this.measuredFrequencies = new ArrayList<Double>();
		this.assigned = false;
	}
	
	/* 
	 * The index mapper is used to find chemical shift in measured peak
	 */
	public Atom (String name, int measuredPeakIndex, int residueNumber, String genericAtomName) {
		this.name = name;
		this.measuredPeakIndex = measuredPeakIndex;
		this.expectedPeaks = new ArrayList<Peak>();
		this.expectedPeaksIndexMapper = new ArrayList<Integer>();
		this.measuredFrequencies = new ArrayList<Double>();
		this.indices = new ArrayList<Integer>();
		this.expectedPeakListAssigned = new ArrayList<Integer>();
		this.assigned = false;
	}
	
	/*
	 * This function adds the given expected peak to the expected peak list of the atom. 
	 */
	public void addExpectedPeak(Peak pk, int j){
		pk.setAssigned(false);
		this.expectedPeaks.add(pk);
		this.expectedPeaksIndexMapper.add(j);
	}
	
	public void increasecf2ByOne(){
		this.setXyz(this.getXyz() + 1);
	}
	
	/*
	 * This function adds the given expected peak to the atom.expectedpeaks with the given assignment index.
	 */
	public void addExpectedPeakWithAssignmentIndex(Peak pk, int j, int startIndex){

		pk.setAssignmentIndex(startIndex+pk.getIndex());
		this.expectedPeaks.add(pk);
		this.expectedPeaksIndexMapper.add(j);
	
	}
	
	public void printiEP(){
		System.out.print("iEP of "+this.name+": ");
		for (int i = 0; i< expectedPeakListAssigned.size(); i++)
			System.out.print(this.expectedPeakListAssigned.get(i) + " ");
		System.out.println();
	}
	
	
	public int getIndicesSingleSpectrum(int spec){
		return this.indices.get(spec);
	}
	
	
	public int getExpectedPeakSize(){
		return this.expectedPeaks.size();
	}
	
	public ArrayList<Peak> getExpectedPeaks(){
		return this.expectedPeaks;
	}
	
	public Peak getExpectedPeak(int i){
		return this.expectedPeaks.get(i);
	}
	
	public int getMeasuredPeakIndices(){
		return this.measuredPeakIndex;
	}
	
	
	public void setMeanDev (double value) {
		this.deviationMeanValue = value;
	}
	
	public double getMeanDev() {
		return this.deviationMeanValue;
	}
	

	/*
	 * This function resets the measures peak frequency values assigned to the atom.
	 * 
	 */
	public void resetMeasuredFrequencyValu(){
		this.assignedMeasuredPeakCount = 0;  
		this.measuredPeakFrequencySum = 0; 
		this.measuredPeakAverage = 0; 
	}
	
	/*
	 * This function resets the measured frequency counters and average for the atom. 
	 */
	public void resetMeasFSpec(){
		
		this.assignedMeasuredPeakCountSpectrum = 0; 
		this.measuredPeakFrequencySumSpectrum = 0; 
		this.measuredPeakAverageSpectrum = 0;  
		
	}
	
	/*
	 * This function sets the "assigned" value of the expected peak "number" to the "assigned" value
	 * for ex: (2,true) will find the expected peak number 2 of the atom and will set it to true
	 */
	public void setExpectedPeakAssigned(int number, boolean assigned){
		int i, size;
		int peakIndex = -1;
		size = this.getExpectedPeakSize();
		
		try{
			for (i = 0; i < size; i++){
				if (this.expectedPeaks.get(i).getIndex() == number){
					peakIndex = i;
					break;
				}
			}
			if (peakIndex == -1){
				System.out.println("Error in the setExpectedPeakAssigned of Atom.class. The given index of the expected peak does not belong to the atom!");
			} else{
				this.expectedPeaks.get(peakIndex).setAssigned(assigned);
			}
		} finally{
			
		}
		
	}
	
	public int getMeasuredFrequencySize(){
		return this.measuredFrequencies.size();
	}
	
	public Double getMeasuredFrequency(int i){
		return this.measuredFrequencies.get(i);
	}
	
	public void addMeasuredFrequency(double currentFrequency){
		this.measuredFrequencies.add(currentFrequency);
	}
	
	public void clearMeasuredFrequency(){
		this.measuredFrequencies.clear();
	}
	
	public void printAtom(){
		int i;
		System.out.println("---------Print Atom-----------");
		System.out.print(fileNumber+"\t"+name+"\t"+type+"\t");
		for (i = 0; i < coordinates.length; i++)
			System.out.print(coordinates[i]+"\t");
		for (i = 0; i < covalentConnectivityAtoms.length; i++)
			System.out.print(covalentConnectivityAtoms[i]+"\t");
		System.out.println(pseudoAtomNumber);
		System.out.println ("Expected Peaks: ");
		for (i = 0; i <this.expectedPeaks.size(); i++)
			this.expectedPeaks.get(i).printPeakLine();
		System.out.print("Measured Frequencies: ");
		for (i = 0; i < this.measuredFrequencies.size(); i++)
			System.out.print(this.measuredFrequencies.get(i)+" ");
		System.out.println ();
		System.out.println("averageFrequency: "+averageFrequency);

		System.out.println("assigned: " + this.assigned);
		System.out.println("accuracy: "+ accuracy);
		System.out.println("partDev: "+ deviationValue);
		System.out.println("referenceFrequency: "+ referenceFrequency);
		System.out.println("partFreq: "+ frequencyValue);
		System.out.println("index: "+ this.index+" (index in Spectrum.atoms)");
		System.out.println("atomsIndex: "+this.atomsIndex+" (index in chemicalShiftAssignment.atoms)");
		System.out.print ("indices: ");
		for (i = 0; i<this.indices.size(); i++)
			System.out.print(this.indices.get(i) + " ");
		System.out.println("(indices of this atom in each spectrum)");
		
		System.out.println ("Search Range up and low: "+this.lowSearchBoundary +" "+this.highSearchBoundary);
		
		System.out.println("--------------------");
		
	}
	
	public void printMeasuredFrequencies() {
		int i;
		System.out.print("Measured Frequencies of atom "+ this.name+": ");
		for (i = 0; i < this.measuredFrequencies.size(); i++)
			System.out.print(this.measuredFrequencies.get(i)+" ");
		System.out.println();
	}
	
	public String getName(){
		return this.name;
	}
	
	public void setName (String name){
		this.name = name;
	}
	
	public boolean getAssignedOld(){
		return this.assigned;
	}
	
	public boolean getAssigned(int[] assignment) {
		boolean assigned = false;
		int i, assignmentIndex;
		
		for (i = 0; i < this.expectedPeaks.size(); i++) {
			assignmentIndex = this.getExpectedPeak(i).getAssignmentIndex();
			if (assignment[assignmentIndex] > -1) {
				assigned = true;
				return assigned;
			}
		}
		return assigned;
	} 
	
	public void setAssigned(boolean assigned){
		this.assigned = assigned;
	}
	
	public int getIndex (){
		return this.index;
	}
	
	public int getiEPSize (){
		return this.expectedPeakListAssigned.size();
	}
	
	public int getiEPValue(int index){
		
		if (this.expectedPeakListAssigned.size()> index)
			return this.expectedPeakListAssigned.get(index);
		else{
			return -1;
		}
	}
	public int dequeueFromiEP(){
		return this.expectedPeakListAssigned.remove(0);
	}
	
	public void initializeiEP (ArrayList<Peak> expPeaks){
		int i;
		
		for (i =0; i< expPeaks.size(); i++){
			if (newValueForIEP(expPeaks.get(i).getAssignmentIndex()) == true){
				this.expectedPeakListAssigned.add(expPeaks.get(i).getAssignmentIndex());
			} else{			//	not a new value, so do not add
			}
		}
		
	}
	
	private boolean newValueForIEP(int expId) {
		int i;
		
		for (i = 0 ; i < this.expectedPeakListAssigned.size(); i++) {
			if (expId == this.expectedPeakListAssigned.get(i)){
				return false;
			}		
		}
		return true;
	}

	public void setIndex(int index){
		this.index = index;
	}
	

	public int getAtomsIndex (){
		return this.atomsIndex;
	}
	
	public void setAtomsIndex(int atomsIndex){
		this.atomsIndex = atomsIndex;
	}

	public int getResidueNumber(){
		return this.residueNumber;
	}
	
	public void setResidueNumber(int residueNumber){
		this.residueNumber = residueNumber;
	}

	public double getScoreLocal(){
		return this.localScoreValue;
	}

	public void setLocalScore(double d){
		this.localScoreValue = d;
	}
	public String getGenericAtomName(){
		return this.genericAtomName;
	}
	
	public void setGenericAtomName(String genericAtomName){
		this.genericAtomName = genericAtomName;
	}
	
	/*
	 * This function adds given frequency to the atom variables
	 */
	public void addFrequencyToChemShiftAtom(double freq){
		
		this.assignedMeasuredPeakCount += 1;
		this.measuredPeakFrequencySum += freq;
		this.measuredPeakAverage = this.measuredPeakFrequencySum / this.assignedMeasuredPeakCount;
	}

	/*
	 * This function adds given frequency to the atom variables from Spectrum
	 */
	public void addFrequencyToSpectrumAtom(double freq){
		

		this.assignedMeasuredPeakCountSpectrum += 1;
		this.measuredPeakFrequencySumSpectrum += freq;
		this.measuredPeakAverageSpectrum = this.measuredPeakFrequencySumSpectrum / this.assignedMeasuredPeakCountSpectrum;
		
	}

	/*
	 * This function calculates the average frequency of the atom and returns this value.
	 */
	public void calculateAverageFrequency() {
		int i, size;
		double sum = 0;
		
		size = this.measuredFrequencies.size();
		
		for (i = 0; i < size; i++){
			sum += this.measuredFrequencies.get(i);
		}
		this.averageFrequency = sum/size;
		
	}
	
	/*
	 * This function returns the predefined average frequency of the atom.
	 */
	public double getAverageFrequency(){
		return this.averageFrequency;
	}
	
	/*
	 * This function sets the average freqeuncy to the given input value i.
	 */
	public void setAverageFrequency(double i) {
		this.averageFrequency = i;
	}

	
	/*
	 * This function assigns the accuracy values of the atoms. 
	 */
	public void setAccuracyValues() {
		String c;
		c = this.name.substring(0, 1);
		switch (c) {
	        case "N": this.accuracy = 0.40000000000000002;
	        	break;
	        case "H": this.accuracy = 2.9999999999999999E-002;
	        	break;
	        case "C":   this.accuracy = 0.40000000000000002;
	        	break;
	        default:  this.accuracy = -1;
	       		break;
		}
	}
	
	public double getAccuracy(){
		return this.accuracy;
	}
	
	public double getDeviationValue(){
		return this.deviationValue;
	}
	
	public void setDeviationValue(double d){
		this.deviationValue = d;
	}
	
	public void increaseDeviationValue(double d){
		this.deviationValue += d;
	}

	public double getReferenceFrequency() {
		return this.referenceFrequency;
	}

	public void setReferenceFrequency(double referenceFrequency) {
		this.referenceFrequency = referenceFrequency;
	}

	public void setCSTableFreq(Double dev) {
		this.csTableFrequency = dev;
	}
	
	public int getAssignedMeasuredPeakCount(){
		return this.assignedMeasuredPeakCount;
	}
	public void setAssignedMeasuredPeakCount(int i){
		this.assignedMeasuredPeakCount = i;
	}
	public double getMeasPeakFreqSum(){
		return this.measuredPeakFrequencySum;
	}
	public void setMeasPFreqSum(double d){
		this.measuredPeakFrequencySum = d;
	}
	public double getMeasuredPeakAverage(){
		return this.measuredPeakAverage;
	}
	public int getAssignedPeakCountSpectrum(){
		return this.assignedMeasuredPeakCountSpectrum;
	}
	public void setAssignedPeakCountSpectrum(int i){
		this.assignedMeasuredPeakCountSpectrum = i;
	}
	public double getMeasuredPeakFreqSumSpectrum(){
		return this.measuredPeakFrequencySumSpectrum;
		}
	public void setMeasuredPeakFreqSum(double d){
		this.measuredPeakFrequencySumSpectrum = d;;
		}
	public double getMeasuredPeakAvgSum(){
		return this.measuredPeakAverageSpectrum;
	}
	
	
	public double getFrequencyValue(){
		return this.frequencyValue;
	}
	
	public void setFrequencyValue(double s){
		this.frequencyValue = s;
	}
	
	public double getGlobalScore(){
		return this.globalScore;
	}
	
	public void setGlobalScore(double s){
		this.globalScore = s;
	}
	
	public double getCsTableFrequency() {
		return this.csTableFrequency;
	}

	public void setMeasuredFrequencies(ArrayList<Double> measuredFrequenciesIn) {
		int i;
		this.measuredFrequencies.clear();
		for (i = 0; i < measuredFrequenciesIn.size(); i++)
			this.measuredFrequencies.add(measuredFrequenciesIn.get(i));
	}
	
	public ArrayList<Double> getMeasuredFrequencies(){
		return this.measuredFrequencies;
	}
	
	public ArrayList<Integer> getExpectedPeaksIndexMapper(){
		return this.expectedPeaksIndexMapper;
	}
	
	public int getExpectedPeaksIndexMapperOf(int i){
		return this.expectedPeaksIndexMapper.get(i);
	}

	/*
	 * This function sets the index of the spec. spectrum in indices
	 */
	public void setIndicesSingleSpectrum(int spec, int a) {
		this.indices.set(spec, a);
	}

	/*
	 * This function returns the smallest shift value of this atom on the search range. 
	 */
	public double getLowSearchBoundary(){
		return this.lowSearchBoundary;
	}
	
	/*
	 * This function returns the highest shift value of this atom on the search range. 
	 */
	
	public double getHighSearchBoundary(){
		return this.highSearchBoundary;
	}
	
	/*
	 * This function sets the low and up boundaries of the search range for this atom.  
	 */
	public void setSearchBoundary(double low, double up){
		this.lowSearchBoundary = low;
		this.highSearchBoundary = up;
	}

	public void setLowSearchBoundary(double d){
		this.lowSearchBoundary = d;
	}

	public void setHighSearchBoundary(double d){
		this.highSearchBoundary = d;
		
	}


	public void setMeasuredFrSpectrum(double measFSpec) {

		this.measuredPeakAverageSpectrum = measFSpec;
	}

	
	public void setMeasuredPeakAverage(double measF) {
		this.measuredPeakAverage = measF;
	}

	public double getDynamicExpectedFreqDeviation() {
		return dynamicExpectedFrequencyDeviation;
	}

	public void setDynamicExpectedFrequencyDeviation(double expFreqDevDynamic) {
		this.dynamicExpectedFrequencyDeviation = expFreqDevDynamic;
	}

	public double getDynamicExpectedFrequency() {
		return dynamicExpectedFrequency;
	}

	public void setExpectedFrequencyValue(double expFreqValDynamic) {
		this.dynamicExpectedFrequency = expFreqValDynamic;
	}


	public void setDynamicExpectedFrequency(double expFreqDevDynamic, double expFreqValDynamic) {
		this.dynamicExpectedFrequencyDeviation = expFreqDevDynamic;
		this.dynamicExpectedFrequency = expFreqValDynamic;
	}

	public double getMeasDevDynamic() {
		return dynamicMeasuredFrequencyValue;
	}

	public void setMeasDevDynamic(double measDevDynamic) {
		this.dynamicMeasuredFrequencyValue = measDevDynamic;
	}

	public double getMeasFreqDynamic() {
		return dynamicMeasuredFrequency;
	}

	public void setDynamicMeasFreq(double measFreqDynamic) {
		this.dynamicMeasuredFrequency = measFreqDynamic;
	}
	
	public void setDynamicMeasuredFrequency(double measDevDynamic, double measFreqDynamic) {
		this.dynamicMeasuredFrequencyValue = measDevDynamic;
		this.dynamicMeasuredFrequency = measFreqDynamic;
	}
	
	/*
	 * This function returns the low and high search boundary shift values for this atom.
	 */
	public double[] getSearchSpace() {
		double[] range = new double[2];
		
		range[0] = this.lowSearchBoundary;
		range[1] = this.highSearchBoundary;
				
		return range;
	}

	public void setScored(boolean b) {

		this.scoreAssigned = b;
		
	}

	public void nullifyiEP() {
		this.expectedPeakListAssigned.clear();
		
	}

	public int getPhaseNr() {
		return phaseNr;
	}

	public void setPhaseNr(int phaseNr) {
		this.phaseNr = phaseNr;
	}

	public double getAverageScoreLocal() {
		return averageScoreLocal;
	}

	public void setAverageScoreLocal(double avgScore) {
		this.averageScoreLocal = avgScore;
	}

	public double getDeviationLocal() {
		return deviationLocal;
	}

	public void setDeviationLocal(double standDev) {
		this.deviationLocal = standDev;
	}

	public void setProbabilitySum(double value){
		this.probabilitySum = value;
	}

	public double getProbabilitySum(){
		return this.probabilitySum;
	}
	
	public void increaseProbabilitySum (double value){
		this.probabilitySum += value;
	}

	public void initializeScalingFactor(int specCount, int nrSpecTypes) {
		int i;
		this.probabilitySumSpectrum = new double[specCount];
		this.scalingFactorSpectrum = new double[nrSpecTypes];
		for (i = 0; i < specCount; i++)
			this.probabilitySumSpectrum[i] = 0;
		for (i = 0; i < nrSpecTypes; i++)
			this.scalingFactorSpectrum[i] = 0;
		
	}

	public void increaseExpInSpec(int specIndex, double value){
		this.probabilitySumSpectrum[specIndex] += value;
	}
	
	public double getProbabilitySumSpectrum(int specIndex){
		return this.probabilitySumSpectrum[specIndex];
	}
	public void setProbabilitySumSpectrum(int specIndex, double value){
		this.probabilitySumSpectrum[specIndex] = value;
	}

	public void setScalingFactorSpectrum(double[] relSpecType, double sumFacExist) {
		int i;
		
		for (i = 0; i < this.scalingFactorSpectrum.length; i++)
			this.scalingFactorSpectrum[i] = relSpecType[i] / sumFacExist;
		
	}

	public double[] getScalingFactorSpectrum() {
		return this.scalingFactorSpectrum;
	}
	
	public double getScalingFactorSpectrum(int specId) {
		return this.scalingFactorSpectrum[specId];
	}

	public void printExpectedPeaks() {
		int i;
		System.out.print("Expected peaks of atom "+ this.name + ": ");
		for (i = 0; i < this.expectedPeaks.size(); i++){
			System.out.print(this.expectedPeaks.get(i).getAssignmentIndex()+ " ");
		}
		System.out.println();
		
		
	}

	public double getMaxScore() {
		return maximumScore;
	}

	public void setMaxScore(double maxScore) {
		this.maximumScore = maxScore;
	}

	public void clearMeasFValues() {
		// TODO Auto-generated method stub
		
	}

	public int getXyz() {
		return xyz;
	}

	public void setXyz(int value) {
		this.xyz = value;
	}


}



