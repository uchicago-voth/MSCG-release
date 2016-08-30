//
//  interaction_hashing.h
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#ifndef _interaction_hashing_h
#define _interaction_hashing_h

#include <vector>

// Return total possible pair interaction types given a number of CG types
inline unsigned calc_n_distinct_pairs(const int n) {return (n * (n + 1) / 2);} 

// Return total possible three body interaction types given a number of CG types
inline unsigned calc_n_distinct_triples(const int n) {return (n * n * (n + 1) / 2);} 

// Return total possible four body interaction types given a number of CG types
inline unsigned calc_n_distinct_quadruples(const int n) {return (n * n * (n * n + 1)/2);}

// Calculate a hash number from a vector of types involved.
int calc_interaction_hash(const std::vector<int> &types, const int n_cg_types);

// Invert a hash number into a vector of types involved.
// Infers the number of types from the size of 'types'.
void invert_interaction_hash(const int m, const int n_cg_types, std::vector<int> &types);

// Calculate a two-body interaction hash number from the two types of site involved.
int calc_two_body_interaction_hash(int i, int j, const int n_cg_types);
// Invert a two-body interaction hash number of the two types of site involved.
void invert_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j);

// Calculate a three-body interaction hash number from the three types of site involved.
int calc_three_body_interaction_hash(int i, int j, int k, const int n_cg_types);
// Invert a three-body interaction hash number of the three types of site involved.
void invert_three_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k);

// Calculate a four-body interaction hash number from the four types of site involved.
int calc_four_body_interaction_hash(int i, int j, int k, int l, const int n_cg_types);
// Invert a four-body interaction hash number of the four types of site involved.
void invert_four_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k, int &l);

// Search a sorted integer table to obtain the index of a desired value.
int SearchIntTable(const std::vector<unsigned> &a, const unsigned m);

// Find the number of active, actual interactions from an array of activation flags.
int calc_n_active_interactions(const int* const activation_flags_by_type, const int tol_n);

// Set up a cached array of interaction type hashes for fast lookup
void set_up_interaction_type_hash_array(int* const activation_flags_by_type, const int n_possible_interactions, std::vector<unsigned> &hash_array);

#endif
