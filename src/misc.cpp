//
//  misc.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cstdlib>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <sstream>

#include "misc.h"

// Global variable assignments

const double VERYLARGE = 1000.0;
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
        fprintf(stderr, "Failed to open file %s.\n", file_name);
        fflush(stderr);
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

void make_negative(std::vector<double> &force_vals)
{
  for(int k = force_vals.size() - 1; k >= 0; k--)
    {
      force_vals[k] = -force_vals[k];
    }
}


// Function to pad 2 vectors so that the first runs between low and high values with fpad

void pad_values_front(const double low, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double fpad)
{
	std::vector<double>::iterator axis_it;
  	std::vector<double>::iterator force_it;
 
	double spacing = axis_vals[1] - axis_vals[0];
	while (axis_vals[0] - spacing > low) {
		axis_it = axis_vals.begin();
 		force_it = force_vals.begin();
 
		axis_vals.insert(axis_it, axis_vals[0] - spacing);
		force_vals.insert(force_it, fpad);	
	}
	
	if (axis_vals[0] - VERYSMALL_F > low) {
		axis_it = axis_vals.begin();
 		force_it = force_vals.begin();
 
		axis_vals.insert(axis_it, low);
		force_vals.insert(force_it, fpad);
	}
}

void pad_values_back(const double high, std::vector<double>& axis_vals, std::vector<double>& force_vals, const double fpad)
{
  double spacing = axis_vals[2] - axis_vals[1];
  int last = axis_vals.size() - 1;
	
  while(force_vals.size() > axis_vals.size()) {
  	force_vals.pop_back();
  }

  while (axis_vals[last] + spacing < high) {
	axis_vals.push_back(axis_vals[last] + spacing);
	force_vals.push_back(fpad);
	last++;
  }
	
  if (axis_vals[last] + VERYSMALL_F < high) { 
	axis_vals.push_back(high);
	force_vals.push_back(fpad);
  }	
}

// Function to pad vectors so that the first runs between low and high values after checking
// that the signs and slopes are reasonable for the front and back.

int pad_values_front_with_fix(std::vector<double>& axis_vals, std::vector<double>& force_vals)
{
  std::vector<double>::iterator axis_it;
  std::vector<double>::iterator force_it;
  int last = axis_vals.size() - 1;
  int flag = 1;
  double spacing = axis_vals[1] - axis_vals[0];
  int i = 0;

  // Find a positive value
  while(force_vals[i] < 0)
    {
      i++;
      if (i  >= last) {
      	// We are out of room
      	flag = -1;
      	i = last - 1;
      	break;
      }
    }
    
  // Keep going until it is non-decreasing
  while(force_vals[i]<force_vals[i+1])
    {
      i++;
      if (i  >= last - 1) {
      	// we are out of room
      	i = last - 1;
      	flag = -1;
      	break;
      }
    }

  // If it failed to meet these requirements (it ran out of room),
  // print a warning to the user 
  // and try again without the positivity requirement.
  if (flag == -1) {
  	i = 0;
  	// Go until it is non-decreasing
  	while(force_vals[i]<force_vals[i+1])
    {
      i++;
      if (i  >= last - 1) {
      	// we are out of room
      	// this time, abort padding front
      	return flag;
      }
    }
  }
  
  // Now, pad interaction
  // first filling in the existing spaces
  while(i > 0)
    {
      force_vals[i-1]=2*force_vals[i]-force_vals[i+1];
      i--;
    }
  // And then adding new ones
  while(axis_vals[0] - spacing > spacing)
    {
      axis_it = axis_vals.begin();
      force_it = force_vals.begin();

      axis_vals.insert(axis_it, axis_vals[0] - spacing);
      force_vals.insert(force_it, 2*force_vals[0] - force_vals[1]);
    }
  return flag;
}

int pad_values_back_with_fix(const double high, std::vector<double>& axis_vals, std::vector<double>& force_vals)
{
  double spacing = axis_vals[2] - axis_vals[1];
  int last = axis_vals.size() - 1;
  int i = last;
  int flag = 1;
 
  if(axis_vals[last] > high) {
  	return flag;
 }
 
  while(force_vals.size() > axis_vals.size()) {
  	force_vals.pop_back();
  }

  // Find a negative value
  while(force_vals[i]>0)
    {
      i--;
      if (i < 1) {
      	i = 1;
      	flag = -1;
      	break;
      }
    }
  // Keep going until it is non-increasing
  while(force_vals[i]>force_vals[i-1])
    {
      i--;
      if (i < 1) {
      	i = 1;
      	flag = -1;
      	break;
      	}
    }

  // If it failed to meet these requirements (it ran out of room),
  // print a warning to the user 
  // and try again without the negativity requirement.
  if (flag == -1) {
  	i = last;
  	while(force_vals[i]>force_vals[i-1]) {
      i--;
      if (i < 1) {
      	// we are out of room
      	// this time, abort padding back
      	return flag;
      }
    }
  }
  
  
  // Now, pad interaction
  // first filling in the existing spaces
  while(i<last)
    {
      force_vals[i+1]= 2.0*force_vals[i] -  force_vals[i-1];
      i++;
    }
  // And then adding new ones
  while (axis_vals[last] + spacing < high)
    {
      axis_vals.push_back(axis_vals[last] + spacing);
      force_vals.push_back(2.0*force_vals[last]  - force_vals[last-1]);
      last++;
    }
  return flag;
}

// Add two sets of "forces" based on their axis values

void add_force_vals(const std::vector<double> &axis_vals, std::vector<double> &force_vals, const std::vector<double> &tab_axis_vals, const std::vector<double> &tab_force_vals)
{
	// Search for the first common axis value between the splines
	int axis_index = 0;
	int tab_index  = 0;
	int last_axis  = axis_vals.size() - 1;
	int last_tab   = tab_axis_vals.size() - 1;
	
	// If tab is lower
	while (axis_vals[axis_index] > tab_axis_vals[tab_index] + VERYSMALL_F) {
		tab_index++;
		if (tab_index >= last_tab) return; // There is no overlap
	}
	// If axis is lower
	while (axis_vals[axis_index] + VERYSMALL_F < tab_axis_vals[tab_index]) {
		axis_index++;
		if (axis_index >= last_axis) return; // There is no overlap
	}
	
	// Now add values until the end
	for ( ; (axis_index < last_axis) && (tab_index < last_tab); axis_index++, tab_index++) {
		force_vals[axis_index] += tab_force_vals[tab_index];
	}
}
    		  
// Remove entries with axis values out of the specified range.

void trim_excess_axis(const double low_value, const double high_value, std::vector<double> &axis_vals, std::vector<double> &force_vals)
{
	int last = axis_vals.size() - 1;
	int first = 0;
	while (axis_vals[last] > VERYSMALL_F + high_value) {
		// remove last entry
		axis_vals.pop_back();
		force_vals.pop_back();
		last--;
	}

	// Check lower value(s)
	while (axis_vals[first] < low_value - VERYSMALL_F) {
		first++;
	}
	if (first != 0) {
		// remove these leading entries
		// only do this once because the operation is rather inefficient
		if (first == 1) {
			//remove a sigle entry
			axis_vals.erase(axis_vals.begin());
			force_vals.erase(force_vals.begin());
		} else {
			// remove a range of entries
			axis_vals.erase(axis_vals.begin(), axis_vals.begin() + first);
			force_vals.erase(force_vals.begin(), force_vals.begin() + first);
		}
	}

}

// Wrap the axis around a boundary if there is more than 1 value to wrap.

void wrap_periodic_axis(const double low_value, const double high_value, std::vector<double> &axis_vals, std::vector<double> &force_vals)
{
	// Determine how many values there are to wrap
	int counter = 0;
	int index = axis_vals.size() - 1;
	while (axis_vals[index] > high_value + VERYSMALL_F) {
		counter++;
		index--;
	}
	
	// Only continue if there is more than 1 value to wrap.
	// If there is only 1 value than it is probably a rounding issue with the spline table.
	// If there are more values, then this interaction range probably passed through a periodic boundary.
	if (counter <= 1) return;
		
	double axis_value, force_value;
	std::vector<double>::iterator axis_it;
  	std::vector<double>::iterator force_it;

	while (counter > 0) {
		// Store last value
		axis_value = axis_vals.back();
		force_value = force_vals.back();
		
		// Remove value
		axis_vals.pop_back();
		force_vals.pop_back();
		
		// Adjust axis
		axis_value -= 360.0;
		
		// Insert value 
		axis_it = axis_vals.begin();
 		force_it = force_vals.begin();
		axis_vals.insert(axis_it, axis_value);
		force_vals.insert(force_it, force_value);	
		
		// adjust counter
		counter--;
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

// A function to match a string against the list of type names.
int match_type(std::string &source, char** name, const int n_types)
{
	char unknown[10];
	sprintf(unknown, "%s", source.c_str());
	// Check if this type is a name.
	for (int i = 0; i < n_types; i++) {
		if( strcmp(unknown, name[i]) == 0 ) {
			return i + 1;
		}
	}
	
	// Check if this type is a number.
	if(isdigit(unknown[0]) != 0) {
		return atoi(unknown); 
	}
	return -1;	
}

void check_and_open_in_stream(std::ifstream &in_stream, const char* filename) 
{
	in_stream.open(filename, std::ifstream::in);
    if (in_stream.fail()) {
		fprintf(stderr, "Problem opening file %s\n", filename);
		fflush(stderr);
		exit(EXIT_FAILURE);
	}
}

void check_and_read_next_line(std::ifstream &in_stream, std::string &line)
{
	if(!std::getline(in_stream, line)) {
			fprintf(stderr, "\nIt appears that the file is no longer open.\n");
			fprintf(stderr, "Please check that you are not attempting to read past the end of the file and try again.\n");
			fflush(stderr);
			exit(EXIT_FAILURE);
	}	
}

void check_and_read_next_line(std::ifstream &in_stream, std::string &line, int &line_num)
{
	check_and_read_next_line(in_stream, line);
	line_num++;	
}

void swap_pair(int& a, int& b) 
{
	int tn = b;
	b = a;
	a = tn;
}
