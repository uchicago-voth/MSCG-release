//
//  types.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdlib>
#include <cmath>
#include <cstdio>

#include "misc.h"

// Global variable assignments

const double VERYSMALL = 1.0e-14;
const float VERYSMALL_F = 1.0e-6; // Small number for single precision
const double MAX_INPUT_FORCE_VALUE = 1000.0; // Filter some noisy data
const int MAX_CG_TYPE_NAME_LENGTH = 24; // Max length for CG type names
const double DEGREES_PER_RADIAN = 180.0 / M_PI;

// An error-catching wrapper for fopen.

FILE* open_file(const char* file_name, const char* mode)
{
    FILE* filepointer = fopen(file_name, mode);
    if (filepointer == NULL) {
        printf("Failed to open file %s.\n", file_name);
        exit(EXIT_FAILURE);
    }
    return filepointer;
}

// Integrate function to calculate a potential from a force and distance vectors.

void integrate_force(const std::vector<double> &axis_vals, const std::vector<double> &force_vals, std::vector<double> &potential_vals) 
{
    potential_vals = std::vector<double>(axis_vals.size());
    potential_vals[axis_vals.size() - 1] = 0.0;
    for(int k = axis_vals.size() - 2; k >= 0; k--) {     
        potential_vals[k] = potential_vals[k+1] + 0.5 * (axis_vals[k + 1] - axis_vals[k]) * (force_vals[k] + force_vals[k + 1]);
    }
}

// Find the index of the minimum value in a vector.

unsigned get_min_index(const std::vector<double> &potential_vals) 
{
    double min_val = potential_vals[0];
    unsigned min_index = 0;
    for(unsigned k = 0; k < potential_vals.size(); k++) {
        if (potential_vals[k] < min_val) {
            min_val = potential_vals[k];
            min_index = k;
        }
    }
    return min_index;
}

// Subtract off minimum value from a vector.

void standardize_potential(std::vector<double> &potential_vals) 
{
    unsigned min_index = get_min_index(potential_vals);
    // Standardize by the minimum value.
    double min_val = potential_vals[min_index];
    for(unsigned k = 0; k < potential_vals.size(); k++) {
        potential_vals[k] -= min_val;
    }
    // run twice to refine any floating-point problems.
    min_val = potential_vals[min_index];
    for(unsigned k = 0; k < potential_vals.size(); k++) {
        potential_vals[k] -= min_val;
    }
}

// C++ tokenizing function for strings.

int StringSplit(std::string source, const char *const delimiter, std::string* results)
{
	int count = 0;
    size_t prev = 0;
    size_t next = 0;

    while ((next = source.find_first_of(delimiter, prev)) != std::string::npos)
    {
        if (next - prev != 0)
        {
            results[count] = source.substr(prev, next - prev);
	        count++;
        }
        prev = next + 1;
    }

    if (prev < source.size())
    {
        results[count] = source.substr(prev);
        count++;
    }
    return count;
}

