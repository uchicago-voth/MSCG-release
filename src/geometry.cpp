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
void subtract_min_image_vectors(const std::array<double, DIMENSION> const &particle_position1, const std::array<double, DIMENSION> const &particle_position2, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement);
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

void subtract_min_image_vectors(const std::array<double, DIMENSION> const &particle_position1, const std::array<double, DIMENSION> const &particle_position2, const real *simulation_box_half_lengths, std::array<double, DIMENSION> &displacement)
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
    std::array<double, DIMENSION> disp02, disp23, disp13;
    int particle_ids_02[2] = {particle_ids[0], particle_ids[2]};
    int particle_ids_23[2] = {particle_ids[2], particle_ids[3]};
    int particle_ids_13[2] = {particle_ids[1], particle_ids[3]};
    subtract_min_image_vectors(particle_ids_02, particle_positions, simulation_box_half_lengths, disp02);
    subtract_min_image_vectors(particle_ids_23, particle_positions, simulation_box_half_lengths, disp23);
    subtract_min_image_vectors(particle_ids_13, particle_positions, simulation_box_half_lengths, disp13);

    // Calculate the angle, which requires many intermediates.
    double rrbc = 1.0 / sqrt(dot_product(disp23, disp23));	// central bond
    std::array<double, DIMENSION> pb, pc, cross_bc;
    cross_product(disp02, disp23, pb);
    cross_product(disp13, disp23, pc);
    cross_product(pb, pc, cross_bc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    //double s = dot_product( disp02, cross_bc) * rpb1 * rpc1 * rrbc; // LAMMPS has a different s
	double s = - dot_product( pb, disp13) * rpb1 * rrbc; // This is the s calculation that LAMMPS used.
    check_sine(s);
    param_val = atan2(s, c) * DEGREES_PER_RADIAN;
    
    // Calculate the derivatives
    double dot02_23 = dot_product(disp02, disp23);
    double dot13_23 = dot_product(disp13, disp23);
    
	for (unsigned i = 0; i < DIMENSION; i++) {
		derivatives[0][i] = - pb[i] * rrbc / pb2;
		derivatives[1][i] =   pc[i] * rrbc / pc2;
		derivatives[2][i] =  (pb[i] * dot02_23 * rrbc / pb2) \
						   - (pc[i] * dot13_23 * rrbc / pc2) \
						   - derivatives[0][i];
	}
    return true;
}

void calc_radius_of_gyration_and_derivatives(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int n_ids, double &param_val, std::array<double, DIMENSION>* &derivatives)
{
	double rg2 = 0.0;
	double rr2;
	std::array<double, DIMENSION> com;
	std::array<double, DIMENSION> displacement;
	int sub_particle_ids[2];
	int id;
	
	// Calculate the center of mass.
	// All positions are wrapped relative to particle_id.
	// The displacements are averaged before adding the offset of the first particle. 
	sub_particle_ids[0] = particle_ids[0];
	for (int j = 0; j < DIMENSION; j++) com[j] = 0.0;
//	printf("REF: %lf, %lf, %lf\n", particle_positions[particle_ids[0]][0], particle_positions[particle_ids[0]][1], particle_positions[particle_ids[0]][2]);
//	printf("COM calculation: \n");
	for (int i = 1; i < n_ids; i++) {
//		printf("\t%d: %lf, %lf, %lf => ", i, particle_positions[particle_ids[i]][0], particle_positions[particle_ids[i]][1], particle_positions[particle_ids[i]][2]);
		sub_particle_ids[1] = particle_ids[i];
		subtract_min_image_vectors(sub_particle_ids, particle_positions, simulation_box_half_lengths, displacement); 
		for (int j = 0; j < DIMENSION; j++) com[j] += displacement[j];
//		printf("\t%d: %lf, %lf, %lf => ",  i, displacement[0],  displacement[1],  displacement[2]);
	}
	for (int i = 0; i < DIMENSION; i++) {
		com[i] /= (double)(n_ids);
		com[i] += particle_positions[particle_ids[0]][i];
	}
//	printf("COM: %lf, %lf, %lf\n\n", com[0], com[1], com[2]);
		
	// Now, simultaneously calculate the derivative and the radius of gyration.
	// The last index is treated separately since its derivative is not explicitly calculated.
	for (int i = 0; i < n_ids - 1; i++) {
		rr2 = 0;
		id = particle_ids[i];
		subtract_min_image_vectors(com, particle_positions[id], simulation_box_half_lengths, displacement);
		for (int j = 0; j < DIMENSION; j++) rr2 += (displacement[j] * displacement[j]);
		for (int j = 0; j < DIMENSION; j++) derivatives[i][j] = 2.0 * displacement[j] / (double)(n_ids);
		rg2 += (rr2 / (double)(n_ids));
	}
	
	// The last index should not have derivatives calculated
	rr2 = 0;
	id = particle_ids[n_ids - 1];
	subtract_min_image_vectors(com, particle_positions[id], simulation_box_half_lengths, displacement);
	for (int j = 0; j < DIMENSION; j++) rr2 += (displacement[j] * displacement[j]);
	rg2 += (rr2 / (double)(n_ids));

	// The scaling is included in the loop over n_ids
	// The returned value is rg (not rg2)
	param_val = sqrt(rg2);	
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
    std::array<double, DIMENSION> disp02, disp23, disp13;
    int particle_ids_02[2] = {particle_ids[0], particle_ids[2]};
    int particle_ids_23[2] = {particle_ids[2], particle_ids[3]};
    int particle_ids_13[2] = {particle_ids[1], particle_ids[3]};
    subtract_min_image_vectors(particle_ids_02, particle_positions, simulation_box_half_lengths, disp02);
    subtract_min_image_vectors(particle_ids_23, particle_positions, simulation_box_half_lengths, disp23);
    subtract_min_image_vectors(particle_ids_13, particle_positions, simulation_box_half_lengths, disp13);

    // Calculate the angle, which requires many intermediates.
    double rrbc = 1.0 / sqrt(dot_product(disp23, disp23));	// central bond
    std::array<double, DIMENSION> pb, pc, cross_bc;
    cross_product(disp02, disp23, pb);
    cross_product(disp13, disp23, pc);
    cross_product(pb, pc, cross_bc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    //double s = dot_product( disp02, cross_bc) * rpb1 * rpc1 * rrbc; // LAMMPS has a different s
	double s = - dot_product( pb, disp13) * rpb1 * rrbc; // This is the s calculation that LAMMPS used.
	check_sine(s);
    param_val = atan2(s, c) * DEGREES_PER_RADIAN;
}

void calc_radius_of_gyration(const int* particle_ids, const std::array<double, DIMENSION>* const &particle_positions, const real *simulation_box_half_lengths, const int num_particles, double &param_val)
{
	double rg2 = 0.0;
	double distance;
	int sub_particle_ids[2];
	int id;
	std::array<double, DIMENSION> com, displacement;
	sub_particle_ids[0] = particle_ids[0];
	for (int j = 0; j < DIMENSION; j++) com[j] = 0.0;
	
	// Calculate the center of mass.
	// All positions are wrapped relative to particle_id.
//	printf("REF: %lf, %lf, %lf\n", particle_positions[particle_ids[0]][0], particle_positions[particle_ids[0]][1], particle_positions[particle_ids[0]][2]);
//	printf("COM calculation: \n");

	for (int i = 1; i < num_particles; i++) {
//		printf("\t%d: %lf, %lf, %lf => ", i, particle_positions[particle_ids[i]][0], particle_positions[particle_ids[i]][1], particle_positions[particle_ids[i]][2]);
		sub_particle_ids[1] = particle_ids[i];
		subtract_min_image_vectors(sub_particle_ids, particle_positions, simulation_box_half_lengths, displacement); 
		for (int j = 0; j < DIMENSION; j++) com[j] += displacement[j];
//		printf("%lf, %lf, %lf \n", displacement[0],  displacement[1],  displacement[2]);
	}
		
	// The displacements are averaged before adding the offset of the first particle. 
	for (int i = 0; i < DIMENSION; i++) {
		com[i] /= (double)(num_particles);
		com[i] += particle_positions[particle_ids[0]][i];
	}
//	printf("COM: %lf, %lf, %lf\n\n", com[0], com[1], com[2]);
	
	// Now, calculate the radius of gyration.
//	printf("RG calculation: \n");
	for (int i = 0; i < num_particles; i++) {
		distance = 0;
		id = particle_ids[i];
//		printf("\t%d: %lf, %lf, %lf => ", i, particle_positions[particle_ids[i]][0], particle_positions[particle_ids[i]][1], particle_positions[particle_ids[i]][2]);
		subtract_min_image_vectors(com, particle_positions[id], simulation_box_half_lengths, displacement);
		for (int j = 0; j < DIMENSION; j++) distance += (displacement[j] * displacement[j]);
		rg2 += distance;
//		printf("%lf, %lf, %lf => %lf\n", displacement[0],  displacement[1],  displacement[2], distance);
	}
	
	// Include the scaling.
	double rg = sqrt( rg2 / (double)(num_particles) );	
//	printf( "%lf = sqrt %lf / %lf\n\n", rg, rg2, (double)(num_particles));
	param_val = rg;
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