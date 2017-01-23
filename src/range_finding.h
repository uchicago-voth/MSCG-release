//
//  range_finding.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _range_finding_h
#define _range_finding_h

struct CG_MODEL_DATA;
struct MATRIX_DATA;

// Initialization of storage for the range value arrays and their computation
void initialize_range_finding_temps(CG_MODEL_DATA* const cg);

// Main output function
void write_range_files(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat);

#endif
