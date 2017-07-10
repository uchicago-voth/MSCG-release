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

#include <fstream>
#include <string>
#include "stdio.h"
#include <vector>

//-------------------------------------------------------------
// Global variable definitions
//-------------------------------------------------------------

extern const double VERYLARGE;
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

// An error-catching wrapper for open c++ filestreams for input
void check_and_open_in_stream(std::ifstream &in_stream, const char* filename);

// An (overloaded) error-catching wrapper for std::getline
void check_and_read_next_line(std::ifstream &in_stream, std::string &line);
void check_and_read_next_line(std::ifstream &in_stream, std::string &line, int &line_num);

// Integrate function to calculate a potential from a force and distance vectors.
void integrate_force(const std::vector<double> &axis_vals, const std::vector<double> &force_vals, std::vector<double> &potential_vals);
void make_negative(std::vector<double> &force_vals);

// Function to pad 2 vectors so that the first runs between low and high values with fpad
void pad_values_front(const double low, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double fpad);
void pad_values_back(const double high, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double fpad);
int pad_values_front_with_fix(std::vector<double>& axis_vals, std::vector<double>& force_vals);
int pad_values_back_with_fix(const double high, std::vector<double>& axis_vals, std::vector<double>& force_vals);

// Add two sets of "forces" based on their axis values
void add_force_vals(const std::vector<double> &axis_vals, std::vector<double> &force_vals, const std::vector<double> &tab_axis_vals, const std::vector<double> &tab_force_vals);

// Remove entries with axis values out of the specified range.
void trim_excess_axis(const double low_value, const double high_value, std::vector<double> &axis_vals, std::vector<double> &force_vals);

// Wrap the axis around a boundary if there is more than 1 value to wrap.
double wrap_periodic_axis(const double low_value, const double high_value, std::vector<double> &axis_vals, std::vector<double> &force_vals);

// Adjust the potential for interactions that wrap around a periodic boundary.
void shift_potential_for_periodicity(const std::vector<double> &axis_vals, const std::vector<double> &force_vals, std::vector<double> &corrected_potential_vals, const double wrapped_axis_value, const double lower_limit, const double upper_limit);

// Find the index of the minimum value in a vector.
unsigned get_min_index(const std::vector<double> &potential_vals);

// Subtract off minimum value from a vector.
void standardize_potential(std::vector<double> &potential_vals);

// A tokenizing function for strings.
int StringSplit(std::string source, const char* const delimiter, std::string* results);

// A function that matches an std::string against type names.
int match_type(std::string &source, char** name, const int n_types);

// A simple function for swapping two numbers.
void swap_pair(int& a, int& b);

#endif
