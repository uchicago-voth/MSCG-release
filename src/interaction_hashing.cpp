//
//  interaction_hashing.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include "interaction_hashing.h"
#include "misc.h"
#include <cmath>
#include <cassert>
#include <cstdio>


// Internal function prototypes.

// Calculate a two-body interaction hash number from the two types of site involved.
int calc_asymmetric_two_body_interaction_hash(int i, int j, const int n_cg_types);
// Invert a two-body interaction hash number of the two types of site involved.
void invert_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j);
// Invert a two-body interaction hash number of the two types of site involved.
void invert_asymmetric_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j);

// Invert a three-body interaction hash number of the three types of site involved.
void invert_three_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k);

// Invert a four-body interaction hash number of the four types of site involved.
void invert_four_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k, int &l);

// Invert a five-body interaction hash number of the four types of site involved.
void invert_five_body_interaction_hash(const int m, const int n_cg_types, int &h, int &i, int &j, int &k, int &l);

int four_body_ij_hash(int i, int j, const int n_cg_types);
int five_body_ij_hash(int i, int j, const int n_cg_types);

// Calculate a hash number from a vector of types involved.
int calc_interaction_hash(const std::vector<int> &types, const int n_cg_types) {
	if (types.size() == 1) {
		return types[0] - 1;
    } else if (types.size() == 2) {
        return calc_two_body_interaction_hash(types[0], types[1], n_cg_types);
    } else if (types.size() == 3) {
        return calc_three_body_interaction_hash(types[0], types[1], types[2], n_cg_types);
    } else if (types.size() == 4) {
        return calc_four_body_interaction_hash(types[0], types[1], types[2], types[3], n_cg_types);
    } else if (types.size() == 5) {
        return calc_five_body_interaction_hash(types[0], types[1], types[2], types[3], types[4], n_cg_types);
    } else {
        assert(false);
        return -1;
    }
}

int calc_asymmetric_interaction_hash(const std::vector<int> &types, const int n_cg_types) {
	if (types.size() == 1) {
		return types[0] - 1;
    } else if (types.size() == 2) {
        return calc_asymmetric_two_body_interaction_hash(types[0], types[1], n_cg_types);
    } else {
    	fprintf(stderr, "asymmetric interaction hash is not designed to handle %d types at the same time!\n", (int)(types.size()));
		assert(false);	
	}
}

// Invert a hash number into a vector of types involved.
// Infers the number of types from the size of 'types'.
void invert_interaction_hash(const int m, const int n_cg_types, std::vector<int> &types) {
    if (types.size() == 1) {
    	types[0] = m + 1;
    } else if (types.size() == 2) {
        invert_two_body_interaction_hash(m, n_cg_types, types[0], types[1]);
    } else if (types.size() == 3) {
        invert_three_body_interaction_hash(m, n_cg_types, types[0], types[1], types[2]);
    } else if (types.size() == 4) {
        invert_four_body_interaction_hash(m, n_cg_types, types[0], types[1], types[2], types[3]);
    } else if (types.size() == 5) {
        invert_five_body_interaction_hash(m, n_cg_types, types[2], types[3], types[4], types[0], types[1]);
    } else {
        assert(false);
    }
}

// Invert a hash number into a vector of types involved.
// Infers the number of types from the size of 'types'.
void invert_asymmetric_interaction_hash(const int m, const int n_cg_types, std::vector<int> &types) {
    if (types.size() == 1) {
    	types[0] = m + 1;
    } else if (types.size() == 2) {
        invert_asymmetric_two_body_interaction_hash(m, n_cg_types, types[0], types[1]);
    } else {
    	fprintf(stderr, "asymmetric interaction hash is not designed to handle %d types at the same time!\n", (int)(types.size()));
		assert(false);	
	}
}

// Calculate a two-body interaction hash number from the two types of site involved.
int calc_two_body_interaction_hash(int i, int j, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types);
    assert(0 < j && j <= n_cg_types);
    if (i > j) swap_pair(i, j);
    return ((i - 1) * n_cg_types + j - 1) - i * (i - 1) / 2;
}

// Calculate a two-body asymmetric interaction hash number from the two types of site involved.
int calc_asymmetric_two_body_interaction_hash(int i, int j, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types);
    assert(0 < j && j <= n_cg_types);
    
    return ((i - 1) * n_cg_types + j - 1);
}

// Invert a two-body interaction hash number to get the two types of site involved.
void invert_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j)
{
    int curr_max_i_hash = calc_two_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m);
    assert(m <= curr_max_i_hash);
    int curr_max_i = n_cg_types;
    while (curr_max_i_hash > m) {
        curr_max_i--;
        curr_max_i_hash -= (n_cg_types - curr_max_i + 1);
    }
    i = curr_max_i;
    j = i + (m - curr_max_i_hash);
}

// Invert a two-body interaction hash number to get the two types of site involved.
void invert_asymmetric_two_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j)
{
    int curr_max_i_hash = calc_two_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m);
    assert(m <= curr_max_i_hash);
    i = m % n_cg_types;
    j = (int)( (m - i) / n_cg_types);
    
    i++;
    j++;
}

// Calculate a three-body interaction hash number from the three types of site involved.
int calc_three_body_interaction_hash(int i, int j, int k, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types); 
	assert(0 < k && k <= n_cg_types);
    if (j > k) swap_pair(j, k);
    return (i - 1) * n_cg_types * (n_cg_types + 1) / 2 + (n_cg_types + 1) * (j - 1) - j * (j + 1) / 2 + k;
}

// Calculate a three-body interaction hash number from the three types of site involved.
void invert_three_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k)
{
	int curr_max_i_hash = calc_three_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m);
    assert(m <= curr_max_i_hash);
    
    i = 1 + floor(m / ((n_cg_types * (n_cg_types + 1)) / 2));
    int remainder_hash = m % ((n_cg_types * (n_cg_types + 1)) / 2);
    invert_two_body_interaction_hash(remainder_hash, n_cg_types, j, k);
}

// Helper function for hashing the first two types of a four type interaction.
int four_body_ij_hash(int i, int j, const int n_cg_types) {
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    int n_ab, n1_ab;
    // The "normal" two body hash with degeneracy correction
    n_ab = ((i - 1) * n_cg_types + j - 1) - i * (i - 1) / 2;
    if (i != j) {
        n1_ab = n_ab - i;
        return n1_ab * n_cg_types * n_cg_types + i * n_cg_types * (n_cg_types + 1) / 2;
    } else {
        n1_ab = n_ab - (i - 1);
        return n1_ab * n_cg_types * n_cg_types + (i - 1) * n_cg_types * (n_cg_types + 1) / 2;
    }
}

// This function is currently not used.
// Helper function for hashing the first two types of a five type interaction.
int five_body_ij_hash(int i, int j, const int n_cg_types) {
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    int n_ab, n1_ab;
    // The "normal" two body hash with degeneracy correction
    n_ab = ((i - 1) * n_cg_types + j - 1) - i * (i - 1) / 2;
    if (i != j) {
        n1_ab = n_ab - i;
        return n1_ab * n_cg_types * n_cg_types * n_cg_types + i * n_cg_types * (n_cg_types + 1) / 2;
    } else {
        n1_ab = n_ab - (i - 1);
        return n1_ab * n_cg_types * n_cg_types * n_cg_types + (i - 1) * n_cg_types * (n_cg_types + 1) / 2;
    }
}

// Calculate a four-body interaction hash number from the four types of site involved.
int calc_four_body_interaction_hash(int i, int j, int k, int l, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    assert(0 < k && k <= n_cg_types); 
	assert(0 < l && l <= n_cg_types);
    int n_ij, n_kl;
    if (i > j) {
    	swap_pair(i, j);
        swap_pair(k, l);
    } else if (i == j) {
    	if (k > l) {
    		swap_pair(k, l);
    	}
    }
    
    n_ij = four_body_ij_hash(i, j, n_cg_types);
    n_kl = (k - 1) * n_cg_types + l - 1;
    if (i == j) {
        if (k > l) {
        	swap_pair(k, l);
        }
        n_kl -= k * (k - 1) / 2;
    }
    return (n_ij + n_kl);
}

void invert_four_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k, int &l)
{
	int curr_max_i_hash = calc_four_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m);
    assert(m <= curr_max_i_hash);
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
        k = int(remainder_hash / n_cg_types) + 1;
        l = (remainder_hash % n_cg_types) + 1;
    } else {
        invert_two_body_interaction_hash(remainder_hash, n_cg_types, k, l);
    }
}


// Calculate a five-body interaction hash number from the four types of site involved.
int calc_five_body_interaction_hash(int i, int j, int k, int l, int m, const int n_cg_types)
{
    assert(0 < i && i <= n_cg_types); 
	assert(0 < j && j <= n_cg_types);
    assert(0 < k && k <= n_cg_types); 
	assert(0 < l && l <= n_cg_types);
	assert(0 < m && m <= n_cg_types);
	
	return calc_two_body_interaction_hash(l, m, n_cg_types);
/*	
    // make sure that the first site has the lower type number
    if (i > m) {
    	swap_pair(i, m);
        swap_pair(j, l);
    }

    int n_ij, n_klm;
    n_ij = four_body_ij_hash(i, j, n_cg_types) * n_cg_types;
    n_klm = calc_three_body_interaction_hash(k, l, m, n_cg_types);
    if (i == j) {
        n_klm += - n_cg_types * k * (k + 1) / 2 - n_cg_types + (k - 1) * n_cg_types * (n_cg_types - 1) / 2 + l + 1; // This does not seem right.
    }
    printf("hash %d\n", n_ij + n_klm);
    return (n_ij + n_klm);
*/
}

void invert_five_body_interaction_hash(const int m, const int n_cg_types, int &i, int &j, int &k, int &l, int &n)
{
	invert_two_body_interaction_hash(m, n_cg_types, l, n);
	i = 0;
	j = 0;
	k = 0;
	/*
	int curr_max_i_hash = calc_five_body_interaction_hash(n_cg_types, n_cg_types, n_cg_types, n_cg_types, n_cg_types, n_cg_types);
    assert(0 <= m);
    assert(m <= curr_max_i_hash);
    int curr_max_i = n_cg_types;
    int curr_max_j = n_cg_types;

    int curr_max_ij_hash = five_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);

    while (curr_max_ij_hash > m) {
        assert(curr_max_j >= curr_max_i);
        if (curr_max_j > curr_max_i) {
            curr_max_j--;
            if (curr_max_j > curr_max_i + 1) {
                // This is the most commonly called (~n_cg_types^2) hash difference & the simplest.
                curr_max_ij_hash -= n_cg_types * n_cg_types;
            } else {
                // Suggestion for improvement: derive the correct simple difference for this update.
                curr_max_ij_hash = five_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);
            }  
        } else {
            curr_max_i--;
            curr_max_j = n_cg_types;
            // Suggestion for improvement: derive the correct simple difference for this update.
            curr_max_ij_hash = five_body_ij_hash(curr_max_i, curr_max_j, n_cg_types);
        }
    }
    i = curr_max_i;
    j = curr_max_j;
    int remainder_hash = m - curr_max_ij_hash;
    invert_three_body_interaction_hash(remainder_hash, n_cg_types, k, l, n);
    */
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

int calc_n_active_interactions(const int* const activation_flags_by_type, const int total_num)
{
    int total_active = 0;
    for (int i = 0; i < total_num; i++) if (activation_flags_by_type[i] == 1) total_active++;
    return total_active;
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