/**
 * @file parsemr.cpp
 *
 * Parser for the open MRI sequence format
 * ---------------------------------------
 *
 * This simple program loads an MR sequence and displays the number of RF, gradient
 * and ADC events. The sequence is defined according to the specification of
 * the open file format available in the docs directory.
 *
 * @author Kelvin Layton <kelvin.layton@uniklinik-freiburg.de>
 */

/** \example parsemr.cpp
 *
 * @brief This simple console program demonstrates how to use the ExternalSequence class.
 *
 * The program performs the following steps:
 * - Constructing an ExternalSequence object.
 * - Setting a custom output function using ExternalSequence::SetPrintFunction such that a # symbol is printed before all messages.
 * - Load a sequence with the ExternalSequence::load member function.
 * - Iterate over sequence blocks and count the occurrence of different events.
 * - Print a summary of the sequence to the standard output.
 *
 */

#include "ExternalSequence.h"

#include <iostream>
#include <fstream>
#include <string>
#include <iomanip>

void custom_print(const std::string &str) {
	std::cout << "# " << str << std::endl;
}

/**
 * @brief Entry point for console program
 */
int main(int argc, char* argv[])
{
	std::string path("");
	if (argc>1) {
		path = std::string(argv[1]);
	}

	// Create sequence object and assign output function
	ExternalSequence seq;
	ExternalSequence::SetPrintFunction(&custom_print);

	// Load sequence file
	if (!seq.load(path)) {
		std::cout << "*** ERROR Cannot load external sequence" << std::endl;
		return 1;
	}

	// Loop through blocks and count events
	int numRf=0, numGx=0, numGy=0, numGz=0, numAdc=0, numDelay=0;
	for (int i=0; i<seq.GetNumberOfBlocks(); i++) {
		SeqBlock* block = seq.GetBlock(i);
		if (block->isADC())   numAdc++;
		if (block->isRF())    numRf++;
		if (block->isDelay()) numDelay++;
		if (block->isTrapGradient(0) || block->isArbitraryGradient(0) ) numGx++;
		if (block->isTrapGradient(1) || block->isArbitraryGradient(1) ) numGy++;
		if (block->isTrapGradient(2) || block->isArbitraryGradient(2) ) numGz++;
		delete block;
	}

	// Display summary of sequence events
	std::cout << std::setw(22) << std::left << "Number of blocks: "     << std::setw(4) << std::right << seq.GetNumberOfBlocks() << std::endl;
	std::cout << std::setw(22) << std::left << "Number of RF pulses: "  << std::setw(4) << std::right << numRf << std::endl;
	std::cout << std::setw(22) << std::left << "Number of GX events: "  << std::setw(4) << std::right << numGx << std::endl;
	std::cout << std::setw(22) << std::left << "Number of GY events: "  << std::setw(4) << std::right << numGy << std::endl;
	std::cout << std::setw(22) << std::left << "Number of GZ events: "  << std::setw(4) << std::right << numGz << std::endl;
	std::cout << std::setw(22) << std::left << "Number of readouts: "   << std::setw(4) << std::right << numAdc << std::endl;
	std::cout << std::setw(22) << std::left << "Number of Delays: "     << std::setw(4) << std::right << numDelay << std::endl;
	std::cout << std::endl;

	return 0;
}

