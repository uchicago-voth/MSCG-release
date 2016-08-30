//
//  interaction_hashing.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include "interaction_hashing.h"
#include <cmath>
#include <cassert>
#include <cstdio>

// Function prototypes for functions used internally.

int four_body_ij_hash(int i, int j, const int n_cg_types);

// Calculate a hash number from a vector of types involved.
int calc_interaction_hash(const std::vector<int> &types, const int n_cg_types) {
    if (types.size() == 2) {
        return calc_two_body_interaction_hash(types[0], types[1], n_cg_types);
    } else if (types.size() == 3) {
        return calc_three_body_interaction_hash(types[0], types[1], types[2], n_cg_types);
    } else if (types.size() == 4) {
        return calc_four_body_interaction_hash(types[0], types[1], types[2], types[3], n_cg_types);
    } else {
        assert(false);
        return -1;
    }
}

// Invert a hash number into a vector of types involved.
// Infers the number of types from the size of 'types'.
void invert_interaction_hash(const int m, const int n_cg_types, std::vector<int> &types) {
    if (types.size() == 2) {
        invert_two_body_interaction_hash(m, n_cg_types, types[0], types[1]);
    } else if (types.size() == 3) {
        invert_three_body_interaction_hash(m, n_cg_types, types[0], types[1], types[2]);
    } else if (types.size() == 4) {
        invert_four_body_interaction_hash(m, n_cg_types, types[0], types[1], types[2], types[3]);
    } else {
        assert(false);
    }
}


// Calculate a two-body interaction hash number from the two types of site involved.
int calc_two_body_interaction_hash(int i, int j, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types);
    assert(0 < j && j <= n_cg_types);
    int tn;
    if (i > j) {
        tn = i;
        i = j;
        j = tn;
    }
    return ((i - 1) * n_cg_types + j - 1) - i * (i - 1) / 2;
}

// Invert a two-body interaction hash number to get the two types of site involved.
void invert_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j)
{
    int curr_max_i_hash = calc_two_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m && m <= curr_max_i_hash);
    int curr_max_i = n_cg_types;
    while (curr_max_i_hash > m) {
        curr_max_i--;
        curr_max_i_hash -= (n_cg_types - curr_max_i + 1);
    }
    i = curr_max_i;
    j = i + (m - curr_max_i_hash);
}

// Calculate a three-body interaction hash number from the three types of site involved.
int calc_three_body_interaction_hash(int i, int j, int k, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types); 
	assert(0 < k && k <= n_cg_types);
    int tn;
    if (j > k) {
        tn = j;
        j = k;
        k = tn;
    }
    return (i - 1) * n_cg_types * (n_cg_types + 1) / 2 + (2 * n_cg_types - j + 2) * (j - 1) / 2 + k - j;
}

// Calculate a three-body interaction hash number from the three types of site involved.
void invert_three_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k)
{
    assert(0 <= m && m <= calc_three_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types, n_cg_types));
    
    i = 1 + floor(m / ((n_cg_types * (n_cg_types + 1)) / 2));
    int remainder_hash = m % ((n_cg_types * (n_cg_types + 1)) / 2);
    invert_two_body_interaction_hash(remainder_hash, n_cg_types, j, k);
}

// Helper function for hashing the first two types of a four type interaction.
int four_body_ij_hash(int i, int j, const int n_cg_types) {
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    int n_ab, n1_ab;
    n_ab = ((i - 1) * n_cg_types + j - 1) - i * (i - 1) / 2;
    if (i != j) {
        n1_ab = n_ab - i;
        return n1_ab * n_cg_types * n_cg_types + i * (n_cg_types * (n_cg_types + 1) / 2);
    } else {
        n1_ab = n_ab - (i - 1);
        return n1_ab * n_cg_types * n_cg_types + (i - 1) * (n_cg_types * (n_cg_types + 1) / 2);
    }
}

// Calculate a four-body interaction hash number from the four types of site involved.
int calc_four_body_interaction_hash(int i, int j, int k, int l, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    assert(0 < k && k <= n_cg_types); 
	assert(0 < l && l <= n_cg_types);
    int tmp, n_ij, n_kl;
    if (i > j) {
            tmp = i;
            i = j;
            j = tmp;
            tmp = k;
            k = l;
            l = tmp;
    }
    n_ij = four_body_ij_hash(i, j, n_cg_types);
    if (i != j) {
        n_kl = (k - 1) * n_cg_types + l - 1;
    } else {
        if (k > l) {
            tmp = k;
            k = l;
            l = tmp;
        }
        n_kl = ((k - 1) * n_cg_types + l - 1) - k * (k - 1) / 2;
    }
    return (n_ij + n_kl);
}

void invert_four_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k, int &l)
{
    assert(0 <= m && m <= calc_four_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types, n_cg_types, n_cg_types));
    int curr_max_i = n_cg_types;
    int curr_max_j = n_cg_types;
    int curr_max_ij_hash = four_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);

    while (curr_max_ij_hash > m) {
        assert(curr_max_j >= curr_max_i);
        if (curr_max_j > curr_max_i) {
            curr_max_j--;
            if (curr_max_j > curr_max_i + 1) {
                // This is the most commonly called (~n_cg_types^2) hash difference & the simplest.
                curr_max_ij_hash -= n_cg_types * n_cg_types;
            } else {
                // Suggestion for improvement: derive the correct simple difference for this update.
                curr_max_ij_hash = four_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);
            }
        } else {
            curr_max_i--;
            curr_max_j = n_cg_types;
            // Suggestion for improvement: derive the correct simple difference for this update.
            curr_max_ij_hash = four_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);
        }
    }
    i = curr_max_i;
    j = curr_max_j;
    int remainder_hash = m - curr_max_ij_hash;
    if (curr_max_j != curr_max_i) {
        k = int(remainder_hash / n_cg_types);
        l = remainder_hash % n_cg_types;
    } else {
        invert_two_body_interaction_hash(remainder_hash, n_cg_types, k, l);
    }
}

// Search a monotonic integer table to obtain the index of a desired value.
int SearchIntTable(const std::vector<unsigned> &a, const unsigned m)
{
    int left, right, mid;
    left = 0;
    right = int(a.size()) - 1;
    int return_val = -1;
    while (left <= right) {
        mid = (left + right) / 2;
        if (m == a[mid]) {
            return_val = mid;
            break;
        } else if (m < a[mid]) right = mid - 1;
        else left = mid + 1;
    }
    return return_val;
}

int calc_n_active_interactions(const int* const activation_flags_by_type, const int tol_n)
{
    int tol = 0;
    for (int i = 0; i < tol_n; i++) if (activation_flags_by_type[i] == 1) tol++;
    return tol;
}

// Set up a cached array of interaction hashes for fast lookup
void set_up_interaction_type_hash_array(int* const activation_flags_by_type, const int n_possible_interactions, std::vector<unsigned> &hash_array)
{
    int i, index_among_actual_interactions = 0;
    for (i = 0; i < n_possible_interactions; i++) {
        if (activation_flags_by_type[i] == 1) {
            hash_array[index_among_actual_interactions] = i;
            index_among_actual_interactions++;
        }
    }
}
