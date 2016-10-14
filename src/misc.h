//
//  misc.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _misc_h
#define _misc_h

typedef float real;
typedef real rvec[3];

#include <cstring>
#include <string>
#include "stdio.h"
#include <vector>

//-------------------------------------------------------------
// Global variable definitions
//-------------------------------------------------------------

extern const double VERYSMALL;
extern const float VERYSMALL_F; // Small number for single precision
extern const double MAX_INPUT_FORCE_VALUE; // Filter some noisy data
extern const int MAX_CG_TYPE_NAME_LENGTH; // Max length for CG type names
extern const double DEGREES_PER_RADIAN;

//-------------------------------------------------------------
// Miscellaneous utility functions
//-------------------------------------------------------------

// An error-catching wrapper for fopen.
FILE* open_file(const char* file_name, const char* mode);

// Integrate function to calculate a potential from a force and distance vectors.
void integrate_force(const std::vector<double> &axis_vals, const std::vector<double> &force_vals, std::vector<double> &potential_vals);

// Find the index of the minimum value in a vector.
unsigned get_min_index(const std::vector<double> &potential_vals);

// Subtract off minimum value from a vector.
void standardize_potential(std::vector<double> &potential_vals);

// A tokenizing function for strings
int StringSplit(std::string source, const char* const delimiter, std::string* results);
#endif
