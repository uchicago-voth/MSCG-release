//
//  geometry.cpp
//  
//
//  Copyright (c) 2016 The Voth Group at The University of Chicago. All rights reserved.
//

#include <cmath>
#include "geometry.h"

#ifndef MAXFLOAT
#define MAXFLOAT 1.0E10
#endif

// Function prototypes for internal functions.
void subtract_min_image_vectors(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement);
void subtract_min_image_particles(const std::array<double, DIMENSION> &particle_position1, const std::array<double, DIMENSION> &particle_position2, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement);
void cross_product(const std::array<double, DIMENSION> &a, const std::array<double, DIMENSION> &b, std::array<double, DIMENSION> &c);
double dot_product(const std::array<double, DIMENSION> &a, const std::array<double, DIMENSION> &b);
double dot_product(const double* a, const double* b);
inline void check_sine(double &s);
inline void check_cos(double &cos_theta);

//------------------------------------------------------------
// Small helper functions used internally.
//------------------------------------------------------------

void subtract_min_image_vectors(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement)
{
    for (int i = 0; i < DIMENSION; i++) {
        displacement[i] = particle_positions[particle_ids[1]][i] - particle_positions[particle_ids[0]][i];
        if (displacement[i] > simulation_box_half_lengths[i]) displacement[i] -= 2.0 * simulation_box_half_lengths[i];
        else if (displacement[i] < -simulation_box_half_lengths[i]) displacement[i] += 2.0 * simulation_box_half_lengths[i];
    }
}

void subtract_min_image_particles(const std::array<double, DIMENSION> &particle_position1, const std::array<double, DIMENSION> &particle_position2, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement)
{
    for (int i = 0; i < DIMENSION; i++) {
        displacement[i] = particle_position2[i] - particle_position1[i];
        if (displacement[i] > simulation_box_half_lengths[i]) displacement[i] -= 2.0 * simulation_box_half_lengths[i];
        else if (displacement[i] < -simulation_box_half_lengths[i]) displacement[i] += 2.0 * simulation_box_half_lengths[i];
    }
}

void get_minimum_image(const int l, std::array<double, DIMENSION>* const &x, const real *simulation_box_half_lengths)
{
    for (int i = 0; i < DIMENSION; i++) {
        if (x[l][i] < 0) x[l][i] += 2.0 * simulation_box_half_lengths[i];
        else if (x[l][i] >= 2.0 * simulation_box_half_lengths[i]) x[l][i] -= 2.0 * simulation_box_half_lengths[i];
    }
}

// NOTE: Cross product is only defined for 2^n - 1 dimensions (and only 3 dimensions in the code at the moment).
void cross_product(const std::array<double, DIMENSION> &a, const std::array<double, DIMENSION> &b, std::array<double, DIMENSION> &c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double dot_product(const std::array<double, DIMENSION> &a, const std::array<double, DIMENSION> &b)
{
    double t = 0.0;
    for (int i = 0; i < DIMENSION; i++) t += a[i] * b[i];
    return t;
}

double dot_product(const double* a, const double* b)
{
    double t = 0.0;
    for (int i = 0; i < DIMENSION; i++) t += a[i] * b[i];
    return t;
}

//------------------------------------------------------------
// Calculate translation-invariant geometrical functions of 
// n particle positions and n-1 of their derivatives.
//------------------------------------------------------------

// Calculate a squared distance and one derivative.

bool conditionally_calc_squared_distance_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<double, DIMENSION>* &derivatives)
{
    double rr2 = 0.0;
    std::array<double, DIMENSION> displacement;
    subtract_min_image_vectors(particle_ids, particle_positions, simulation_box_half_lengths, displacement);
    for (int i = 0; i < DIMENSION; i++) {
        rr2 += displacement[i] * displacement[i];
    }
    if (rr2 > cutoff2) {
        return false;
    } else {
        for (int i = 0; i < DIMENSION; i++) {
            derivatives[0][i] = 2.0 * displacement[i];
        }
        param_val = rr2;
        return true;
    }
}

// Calculate a distance and one derivative.

bool conditionally_calc_distance_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<double, DIMENSION>* &derivatives)
{
    bool within_cutoff = conditionally_calc_squared_distance_and_derivatives(particle_ids, particle_positions, simulation_box_half_lengths, cutoff2, param_val, derivatives);

    if (within_cutoff) {
        param_val = sqrt(param_val);
        for (int i = 0; i < DIMENSION; i++) {
            derivatives[0][i] = 0.5 * derivatives[0][i] / param_val;
        }
        return true;
    } else {
        return false;
    }
}

// Calculate the angle between three particles and its derivatives.

bool conditionally_calc_angle_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<double, DIMENSION>* &derivatives)
{   
    std::array<double, DIMENSION>* dist_derivs_20 = new std::array<double, DIMENSION>[1];
    std::array<double, DIMENSION>* dist_derivs_21 = new std::array<double, DIMENSION>[1];
    int particle_ids_20[2] = {particle_ids[2], particle_ids[0]};
    int particle_ids_21[2] = {particle_ids[2], particle_ids[1]};
    double rr2_20, rr2_21;
    bool within_cutoff_20 = conditionally_calc_squared_distance_and_derivatives(particle_ids_20, particle_positions, simulation_box_half_lengths, cutoff2, rr2_20, dist_derivs_20);
    bool within_cutoff_21 = conditionally_calc_squared_distance_and_derivatives(particle_ids_21, particle_positions, simulation_box_half_lengths, cutoff2, rr2_21, dist_derivs_21);
    
    if (!within_cutoff_20 || !within_cutoff_21) {
    	delete [] dist_derivs_20;
    	delete [] dist_derivs_21;
        return false;
    } else {
        // Calculate the cosine
        double rr_20 = sqrt(rr2_20);
        double rr_21 = sqrt(rr2_21);
		double cos_theta = dot_product(dist_derivs_20[0], dist_derivs_21[0]) / (4.0 * rr_20 * rr_21);
		check_cos(cos_theta);
        
        // Calculate the angle.
        double theta = acos(cos_theta);
        param_val = theta * DEGREES_PER_RADIAN;

        // Calculate the derivatives.
        double sin_theta = sin(theta);
        double rr_01_1 = 1.0 / (rr_20 * rr_21 * sin_theta);
        double rr_00c = cos_theta / (rr_20 * rr_20 * sin_theta);
        double rr_11c = cos_theta / (rr_21 * rr_21 * sin_theta);

        for (unsigned i = 0; i < DIMENSION; i++) {
        	// derivatives for the end particles
        	derivatives[0][i] = 0.5 * (dist_derivs_21[0][i] * rr_01_1 - rr_00c * dist_derivs_20[0][i]);
            derivatives[1][i] = 0.5 * (dist_derivs_20[0][i] * rr_01_1 - rr_11c * dist_derivs_21[0][i]);
        }
        delete [] dist_derivs_20;
    	delete [] dist_derivs_21;
        return true;
    }
}

// Calculate a the cosine of an angle along with its derivatives.

bool conditionally_calc_angle_and_intermediates(const int* particle_ids, std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const double cutoff2, std::array<double, DIMENSION>* &dist_derivs_20, std::array<double, DIMENSION>* &dist_derivs_21, std::array<double, DIMENSION>* &derivatives, double &param_val, double &rr_20, double &rr_21)
{
    int particle_ids_20[2] = {particle_ids[2], particle_ids[0]};
    int particle_ids_21[2] = {particle_ids[2], particle_ids[1]};
    double rr2_20, rr2_21;
    bool within_cutoff_20 = conditionally_calc_squared_distance_and_derivatives(particle_ids_20, particle_positions, simulation_box_half_lengths, cutoff2, rr2_20, dist_derivs_20);
    bool within_cutoff_21 = conditionally_calc_squared_distance_and_derivatives(particle_ids_21, particle_positions, simulation_box_half_lengths, cutoff2, rr2_21, dist_derivs_21);
    
    if (!within_cutoff_20 || !within_cutoff_21) {
        return false;
    } else {
        // Calculate the cosine
        rr_20 = sqrt(rr2_20);
        rr_21 = sqrt(rr2_21);
		double cos_theta = dot_product(dist_derivs_20[0], dist_derivs_21[0]) / (4.0 * rr_20 * rr_21);
        check_cos(cos_theta);
        
        // Calculate the angle.
        double theta = acos(cos_theta);
        param_val = theta * DEGREES_PER_RADIAN;

        // Calculate the derivatives.
        double sin_theta = sin(theta);
        double rr_01_1 = 1.0 / (rr_20 * rr_21 * sin_theta);
        double rr_00c = cos_theta / (rr_20 * rr_20 * sin_theta);
        double rr_11c = cos_theta / (rr_21 * rr_21 * sin_theta);

        for (unsigned i = 0; i < DIMENSION; i++) {
            derivatives[0][i] = 0.5 * (dist_derivs_21[0][i] * rr_01_1 - rr_00c * dist_derivs_20[0][i]);
            derivatives[1][i] = 0.5 * (dist_derivs_20[0][i] * rr_01_1 - rr_11c * dist_derivs_21[0][i]);
        }
    }    
    return true;
}

// Calculate a dihedral angle and its derivatives.

bool conditionally_calc_dihedral_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<double, DIMENSION>* &derivatives)
{
    // Find the relevant displacements for defining the angle.
    std::array<double, DIMENSION> disp03, disp23, disp12;
    int particle_ids_03[2] = {particle_ids[3], particle_ids[0]};
    int particle_ids_23[2] = {particle_ids[3], particle_ids[2]};
    int particle_ids_12[2] = {particle_ids[2], particle_ids[1]};
    subtract_min_image_vectors(particle_ids_03, particle_positions, simulation_box_half_lengths, disp03);
    subtract_min_image_vectors(particle_ids_23, particle_positions, simulation_box_half_lengths, disp23);
    subtract_min_image_vectors(particle_ids_12, particle_positions, simulation_box_half_lengths, disp12);

    // Calculate the dihedral, which requires many intermediates.
    // It is the dot product of the cross product. 
    // Note: To calculate the cosine, the vectors in the final dot product need to be effectively normalized.
    double rrbc = 1.0 / sqrt(dot_product(disp23, disp23));	// central bond
    std::array<double, DIMENSION> pb, pc, cross_bc;
    cross_product(disp03, disp23, pb);
    cross_product(disp12, disp23, pc);
    cross_product(pb, pc, cross_bc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    //double s = dot_product( disp03, cross_bc) * rpb1 * rpc1 * rrbc; // LAMMPS has a different s
	double s = - dot_product( pb, disp12) * rpb1 * rrbc; // This is the s calculation that LAMMPS used.
    check_sine(s);
    param_val = atan2(s, c) * DEGREES_PER_RADIAN;
    
    // Calculate the derivatives
    double dot03_23 = dot_product(disp03, disp23);
    double dot12_23 = dot_product(disp12, disp23);
    
	for (unsigned i = 0; i < DIMENSION; i++) {
		derivatives[0][i] = - pb[i] * rrbc / pb2;
		derivatives[1][i] =   pc[i] * rrbc / pc2;
		derivatives[2][i] =  (pb[i] * dot03_23 * rrbc / pb2) \
						   - (pc[i] * dot12_23 * rrbc / pc2) \
						   - derivatives[0][i];
	}
    return true;
}

void calc_radius_of_gyration_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int num_particles, double &param_val, std::array<double, DIMENSION>* &derivatives)
{

	double rg2 = 0.0;
	int sub_particle_ids[2];
	std::array<double, DIMENSION>* whole_molecule = new std::array<double, DIMENSION>[num_particles];
	std::array<double, DIMENSION> com, displacement;
	
	// Get individual particle displacements to reconstruct the whole molecule positions continuously
	for (int j = 0; j < DIMENSION; j++) whole_molecule[0][j] = 0.0;
	for (int i = 1; i < num_particles; i++) {
		sub_particle_ids[0] = particle_ids[i-1];
		sub_particle_ids[1] = particle_ids[i];
		subtract_min_image_vectors(sub_particle_ids, particle_positions, simulation_box_half_lengths, displacement); 
		for (int j = 0; j < DIMENSION; j++) whole_molecule[i][j] = whole_molecule[i - 1][j] + displacement[j];		
	}
	
	// Calculate the center of mass.
	// The displacements are averaged relative to the first particle (which is the center of this reference frame). 
	for (int j = 0; j < DIMENSION; j++) com[j] = 0.0;
	for (int i = 0; i < num_particles; i++) {
		for (int j = 0; j < DIMENSION; j++) com[j] += whole_molecule[i][j];
	}		
	for (int j = 0; j < DIMENSION; j++) com[j] /= (double)(num_particles);

	// Now, simultaneously calculate the derivative and the radius of gyration.
	// The last index is treated separately since its derivative is not explicitly calculated.
	for (int i = 0; i < num_particles - 1; i++) {
		for (int j = 0; j < DIMENSION; j++) displacement[j] = whole_molecule[particle_ids[i]][j] - com[j];
		for (int j = 0; j < DIMENSION; j++) rg2 += (displacement[j] * displacement[j]);
		for (int j = 0; j < DIMENSION; j++) derivatives[i][j] = 2.0 * displacement[j] / (double)(num_particles);
	}

	// The last index should not have derivatives calculated
	for (int j = 0; j < DIMENSION; j++) displacement[j] = whole_molecule[particle_ids[num_particles - 1]][j] - com[j];
	for (int j = 0; j < DIMENSION; j++) rg2 += (displacement[j] * displacement[j]);

	// Include the scaling.
	param_val = rg2 / (double)(num_particles);	
}

void calc_fraction_helical_and_derivatives(const int* particle_ids, std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int n_ids, double &param_val, std::array<double,DIMENSION>* &derivatives, const int* helical_ids, const int n_helical_ids, const int r0, const int sigma2)
{
	std::array<double, DIMENSION>* displacement = new std::array<double, DIMENSION>[1];
	int sub_particle_ids[2];
	double distance, diff, contribution, deriv_magnitude;
	param_val = 0.0;
	
	// zero derivatives to start
	for (int i = 0; i < n_ids - 1; i++) {
		for (int k = 0; k < DIMENSION; k++) {
			derivatives[i][k] = 0.0;
		}
	}
	
	for (int i = 0; i < n_helical_ids; i++) {
		sub_particle_ids[0] = helical_ids[2 * i];
		sub_particle_ids[1] = helical_ids[2 * i + 1];
		conditionally_calc_distance_and_derivatives(sub_particle_ids, particle_positions, simulation_box_half_lengths, 1000000.0, distance, displacement);

		diff = distance - r0;
		contribution = exp( - 0.5 * diff * diff / sigma2);
		param_val += contribution;

		deriv_magnitude = diff * contribution / (sigma2 * distance * n_helical_ids);
		
		// Accumulate this derivative as long as the first index is not the last particle in particle_ids
		if (sub_particle_ids[0] != particle_ids[n_ids - 1]) {
			// find this index in particle_ids
			int index = -1;
			for (int j = 0; j < n_ids - 1; j++) {
				if (particle_ids[j] == sub_particle_ids[0]) {
					index = j;
					break;
				}	
			}
			if (index == -1) {
				printf("Error in indices for helical calculation!\n");
			} else {
				for (int k = 0; k < DIMENSION; k++) {
					derivatives[index][k] += displacement[0][k] * deriv_magnitude;
				}
			}
		}

		// Accumulate this derivative as long as the second index is not the last particle in particle_ids
		if (sub_particle_ids[1] != particle_ids[n_ids - 1]) {
			// find this index in particle_ids
			int index = -1;
			for (int j = 0; j < n_ids - 1; j++) {
				if (particle_ids[j] == sub_particle_ids[1]) {
					index = j;
					break;
				}	
			}
			if (index == -1) {
				printf("Error in indices for helical calculation!\n");
			} else {
				for (int k = 0; k < DIMENSION; k++) {
					derivatives[index][k] -= displacement[0][k] * deriv_magnitude;
				}
			}		
		}
		
	}
	param_val /= (double)(n_helical_ids);
	delete [] displacement;
}

//------------------------------------------------------------
// Without derivatives.
//------------------------------------------------------------

// Calculate a squared distance.

void calc_squared_distance(const int* particle_ids, const std::array<double, DIMENSION>*const &particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
	std::array<double, DIMENSION> displacement;
	double rr2 = 0.0;
    for (int i = 0; i < DIMENSION; i++) {
        subtract_min_image_vectors(particle_ids, particle_positions, simulation_box_half_lengths, displacement);  
        rr2 += displacement[i] * displacement[i];
    }
    param_val = rr2;
}

// Calculate a distance.

void calc_distance(const int* particle_ids, std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
    calc_squared_distance(particle_ids, particle_positions, simulation_box_half_lengths, param_val);
    param_val = sqrt(param_val);
}

// Calculate the angle between three particles.

void calc_angle(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, double &param_val)
{   
    std::array<double, DIMENSION>* dist_derivs_20 = new std::array<double, DIMENSION>[1];
    std::array<double, DIMENSION>* dist_derivs_21 = new std::array<double, DIMENSION>[1];
    int particle_ids_20[2] = {particle_ids[2], particle_ids[0]};
    int particle_ids_21[2] = {particle_ids[2], particle_ids[1]};
    double rr2_20, rr2_21;
    conditionally_calc_squared_distance_and_derivatives(particle_ids_20, particle_positions, simulation_box_half_lengths, MAXFLOAT, rr2_20, dist_derivs_20);
    conditionally_calc_squared_distance_and_derivatives(particle_ids_21, particle_positions, simulation_box_half_lengths, MAXFLOAT, rr2_21, dist_derivs_21);
    
    // Calculate the cosine
    double rr_20 = sqrt(rr2_20);
    double rr_21 = sqrt(rr2_21);
    double cos_theta = dot_product(dist_derivs_20[0], dist_derivs_21[0]) / (4.0 * rr_20 * rr_21);
    check_cos(cos_theta);
    
    // Calculate the angle.
    double theta = acos(cos_theta);
    param_val = theta * DEGREES_PER_RADIAN;

    delete [] dist_derivs_20;
    delete [] dist_derivs_21;
}

// Calculate a dihedral angle.

void calc_dihedral(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
    // Find the relevant displacements for defining the angle.
    std::array<double, DIMENSION> disp03, disp23, disp12;
    int particle_ids_03[2] = {particle_ids[3], particle_ids[0]};
    int particle_ids_23[2] = {particle_ids[3], particle_ids[2]};
    int particle_ids_12[2] = {particle_ids[2], particle_ids[1]};
    subtract_min_image_vectors(particle_ids_03, particle_positions, simulation_box_half_lengths, disp03);
    subtract_min_image_vectors(particle_ids_23, particle_positions, simulation_box_half_lengths, disp23);
    subtract_min_image_vectors(particle_ids_12, particle_positions, simulation_box_half_lengths, disp12);

    // Calculate the dihedral, which requires many intermediates.
    // It is the dot product of the cross product. 
    // Note: To calculate the cosine, the vectors in the final dot product need to be effectively normalized.
    double rrbc = 1.0 / sqrt(dot_product(disp23, disp23));	// central bond
    std::array<double, DIMENSION> pb, pc, cross_bc;
    cross_product(disp03, disp23, pb);
    cross_product(disp12, disp23, pc);
    cross_product(pb, pc, cross_bc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    //double s = dot_product( disp03, cross_bc) * rpb1 * rpc1 * rrbc; // LAMMPS has a different s
	double s = - dot_product( pb, disp12) * rpb1 * rrbc; // This is the s calculation that LAMMPS used.
	check_sine(s);
    param_val = atan2(s, c) * DEGREES_PER_RADIAN;
}

void calc_radius_of_gyration(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int num_particles, double &param_val)
{
	double rg2 = 0.0;
	int sub_particle_ids[2];
	std::array<double, DIMENSION>* whole_molecule = new std::array<double, DIMENSION>[num_particles];
	std::array<double, DIMENSION> com, displacement;
	
	// Get individual particle displacements to reconstruct the whole molecule positions continuously
	for (int j = 0; j < DIMENSION; j++) whole_molecule[0][j] = 0.0;
	for (int i = 1; i < num_particles; i++) {
		sub_particle_ids[0] = particle_ids[i-1];
		sub_particle_ids[1] = particle_ids[i];
		subtract_min_image_vectors(sub_particle_ids, particle_positions, simulation_box_half_lengths, displacement); 
		for (int j = 0; j < DIMENSION; j++) whole_molecule[i][j] = whole_molecule[i - 1][j] + displacement[j];		
	}
	
	// Calculate the center of mass.
	// The displacements are averaged relative to the first particle (which is the center of this reference frame). 
	for (int j = 0; j < DIMENSION; j++) com[j] = 0.0;
	for (int i = 0; i < num_particles; i++) {
		for (int j = 0; j < DIMENSION; j++) com[j] += whole_molecule[i][j];
	}		
	for (int j = 0; j < DIMENSION; j++) com[j] /= (double)(num_particles);
	
	// Now, calculate the radius of gyration as distances from the center of mass.
	for (int i = 0; i < num_particles; i++) {
		for (int j = 0; j < DIMENSION; j++) displacement[j] = whole_molecule[particle_ids[i]][j] - com[j];
		for (int j = 0; j < DIMENSION; j++) rg2 += (displacement[j] * displacement[j]);
	}
	
	// Include the scaling.
	param_val = rg2 / (double)(num_particles);	
}


void calc_fraction_helical(const int* particle_ids, std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int num_particles, double &param_val, const int* helical_ids, const int n_helical_ids, const int r0, const int sigma2)
{
	int sub_particle_ids[2];
	double distance, diff;
	param_val = 0.0;
	
	for (int i = 0; i < n_helical_ids; i++) {
		sub_particle_ids[0] = helical_ids[2 * i];
		sub_particle_ids[1] = helical_ids[2 * i + 1];
		calc_distance(sub_particle_ids, particle_positions, simulation_box_half_lengths, distance);
		diff = distance - r0;
		param_val += exp( - 0.5 * diff * diff / sigma2);
	}
	param_val /= (double)(n_helical_ids);
}

inline void check_sine(double &s)
{
    if (s < 0.0 && s > -VERYSMALL_F) s = -VERYSMALL_F;
    if (s > 0.0 && s < VERYSMALL_F) s = VERYSMALL_F;
}

inline void check_cos(double &cos_theta)
{
        double max = 1.0 - VERYSMALL_F;
        double min = -1.0 + VERYSMALL_F;
        if (cos_theta > max) cos_theta = max;
        else if (cos_theta < min) cos_theta = min;
}