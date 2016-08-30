//
//  batch_fm_combination.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _batch_fm_combination_h
#define _batch_fm_combination_h

struct CG_MODEL_DATA;
struct MATRIX_DATA;

// Specialized deallocation for batch-processing; fewer matrix building
// temporaries are allocated than in comparable functions elsewhere.

void free_AV_fm_setup_temp_variables(CG_MODEL_DATA* const cg, MATRIX_DATA* const mat);

// An extra solver for accumulation matrices using a more sophisticated 
// regularization scheme. Not currently used anywhere.

void solve_regularized_accumulation_form_fm_equations(MATRIX_DATA* const mat);

#endif
