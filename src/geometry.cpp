//
//  geometry.cpp
//  
//
//  Created by James Dama on 12/14/15 from force_computation.h by James Dama
//  Copyright (c) 2015 University of Chicago. All rights reserved.
//

#include <cmath>
#include "geometry.h"

#ifndef MAXFLOAT
#define MAXFLOAT 1.0E10
#endif


// Function prototypes for functions used internally.

void subtract_min_image_vectors(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, std::array<double, 3> &displacement);
void cross_product(const std::array<double, 3> &a, const std::array<double, 3> &b, std::array<double, 3> &c);
double dot_product(const std::array<double, 3> &a, const std::array<double, 3> &b);

//------------------------------------------------------------
// Small helper functions used internally.
//------------------------------------------------------------

void subtract_min_image_vectors(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, std::array<double, 3> &displacement)
{
    for (int i = 0; i < 3; i++) {
        displacement[i] = particle_positions[particle_ids[1]][i] - particle_positions[particle_ids[0]][i];
        if (displacement[i] > simulation_box_half_lengths[i]) displacement[i] -= 2.0 * simulation_box_half_lengths[i];
        else if (displacement[i] < -simulation_box_half_lengths[i]) displacement[i] += 2.0 * simulation_box_half_lengths[i];
    }
}

void cross_product(const std::array<double, 3> &a, const std::array<double, 3> &b, std::array<double, 3> &c)
{
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

double dot_product(const std::array<double, 3> &a, const std::array<double, 3> &b)
{
    double t = 0.0;
    for (int i = 0; i < 3; i++) t += a[i] * b[i];
    return t;
}

//------------------------------------------------------------
// Calculate translation-invariant geometrical functions of 
// n particle positions and n-1 of their derivatives.
//------------------------------------------------------------

// Calculate a squared distance and one derivative.

bool conditionally_calc_squared_distance_and_derivatives(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<std::array<double, 3>, 1> &derivatives)
{
    double rr2 = 0.0;
    std::array<double, 3> displacement;
    subtract_min_image_vectors(particle_ids, particle_positions, simulation_box_half_lengths, displacement);
    for (int i = 0; i < 3; i++) {
        rr2 += displacement[i] * displacement[i];
    }
    if (rr2 > cutoff2) {
        return false;
    } else {
        for (int i = 0; i < 3; i++) {
            derivatives[0][i] = 2.0 * displacement[i];
        }
        param_val = rr2;
        return true;
    }
}

// Calculate a distance and one derivative.

bool conditionally_calc_distance_and_derivatives(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<std::array<double, 3>, 1> &derivatives)
{
    bool within_cutoff = conditionally_calc_squared_distance_and_derivatives(particle_ids, particle_positions, simulation_box_half_lengths, cutoff2, param_val, derivatives);

    if (within_cutoff) {
        param_val = sqrt(param_val);
        for (int i = 0; i < 3; i++) {
            derivatives[0][i] = 0.5 * derivatives[0][i] / param_val;
        }
        return true;
    } else {
        return false;
    }
}

// Calculate the angle between three particles and its derivatives.

bool conditionally_calc_angle_and_derivatives(const std::array<int, 3> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<std::array<double, 3>, 2> &derivatives)
{   
    std::array<std::array<double, 3>, 1> dist_derivs_01;
    std::array<std::array<double, 3>, 1> dist_derivs_02;
    double rr2_01, rr2_02;
    bool within_cutoff_01 = conditionally_calc_squared_distance_and_derivatives(std::array<int, 2>({{particle_ids[0], particle_ids[1]}}), particle_positions, simulation_box_half_lengths, cutoff2, rr2_01, dist_derivs_01);
    bool within_cutoff_02 = conditionally_calc_squared_distance_and_derivatives(std::array<int, 2>({{particle_ids[0], particle_ids[2]}}), particle_positions, simulation_box_half_lengths, cutoff2, rr2_02, dist_derivs_02);
    
    if (!within_cutoff_01 || !within_cutoff_02) {
        return false;
    } else {
        // Calculate the cosine
        double rr_01 = sqrt(rr2_01);
        double rr_02 = sqrt(rr2_02);
		double cos_theta = dot_product(dist_derivs_01[0], dist_derivs_02[0]) / (4 * rr_01 * rr_02);
        double max = 1.0 - VERYSMALL_F;
        double min = -1.0 + VERYSMALL_F;
        if (cos_theta > max) cos_theta = max;
        else if (cos_theta < min) cos_theta = min;
        
        // Calculate the angle.
        double theta = acos(cos_theta);
        param_val = theta * DEGREES_PER_RADIAN;

        // Calculate the derivatives.
        double sin_theta = sin(theta);
        double rr_12_1 = 1.0 / (rr_01 * rr_02 * sin_theta);
        double rr_11c = cos_theta / (rr_01 * rr_01 * sin_theta);
        double rr_22c = cos_theta / (rr_02 * rr_02 * sin_theta);

        for (unsigned i = 0; i < 3; i++) {
            derivatives[0][i] = 0.5 * (-dist_derivs_02[0][i] * rr_12_1 + rr_11c * dist_derivs_01[0][i]);
            derivatives[1][i] = 0.5 * (-dist_derivs_01[0][i] * rr_12_1 + rr_22c * dist_derivs_02[0][i]);
        }
        return true;
    }
}

// Calculate a dihedral angle and its derivatives.

bool conditionally_calc_dihedral_and_derivatives(const std::array<int, 4> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, const double cutoff2, double &param_val, std::array<std::array<double, 3>, 3> &derivatives)
{
    // Find the relevant displacements for defining the angle.
    std::array<double, 3> disp20, disp01, disp23;
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[2], particle_ids[0]}}), particle_positions, simulation_box_half_lengths, disp20);
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[0], particle_ids[1]}}), particle_positions, simulation_box_half_lengths, disp01);
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[2], particle_ids[3]}}), particle_positions, simulation_box_half_lengths, disp23);
    
    // Calculate the angle, which requires many intermediates.
    double rrbc = 1.0 / sqrt(dot_product(disp01, disp01));
    std::array<double, 3> pb, pc;
    cross_product(disp20, disp01, pb);
    cross_product(disp01, disp23, pc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    double s = (disp01[0] * (pc[1] * pb[2] - pc[2] * pb[1]) + disp01[1] * (pb[0] * pc[2] - pb[2] * pc[0]) + disp01[2] * (pc[0] * pb[1] - pc[1] * pb[0])) * (rpb1 * rpc1 * rrbc);
    if (s < 0.0 && s > -VERYSMALL_F) s = -VERYSMALL_F;
    if (s > 0.0 && s < VERYSMALL_F) s = VERYSMALL_F;
    param_val = atan2(s, c);

    // Calculate the derivatives, which requires many more intermediates.
    double gamma = rpb1 * rpc1 / s;
    double cb = pbpc * rpb1 * rpb1;
    double cc = pbpc * rpc1 * rpc1;
    
    std::array<double, 3> w1, w2, w3, w4;
    w1[0] = (-pc[1] * disp01[2] + pc[2] * disp01[1]) - cb * (-pb[1] * disp01[2] + pb[2] * disp01[1]);
    w1[1] = (pc[0] * disp01[2] - pc[2] * disp01[0]) - cb * (pb[0] * disp01[2] - pb[2] * disp01[0]);
    w1[2] = (-pc[0] * disp01[1] + pc[1] * disp01[0]) - cb * (-pb[0] * disp01[1] + pb[1] * disp01[0]);
    
    w2[0] = (-pc[1] * disp20[2] + pc[2] * disp20[1]) - cb * (-pb[1] * disp20[2] + pb[2] * disp20[1]);
    w2[1] = (pc[0] * disp20[2] - pc[2] * disp20[0]) - cb * (pb[0] * disp20[2] - pb[2] * disp20[0]);
    w2[2] = (-pc[0] * disp20[1] + pc[1] * disp20[0]) - cb * (-pb[0] * disp20[1] + pb[1] * disp20[0]);
    
    w3[0] = (-pb[1] * disp23[2] + pb[2] * disp23[1]) - cc * (-pc[1] * disp23[2] + pc[2] * disp23[1]);
    w3[1] = (pb[0] * disp23[2] - pb[2] * disp23[0]) - cc * (pc[0] * disp23[2] - pc[2] * disp23[0]);
    w3[2] = (-pb[0] * disp23[1] + pb[1] * disp23[0]) - cc * (-pc[0] * disp23[1] + pc[1] * disp23[0]);
    
    w4[0] = (-pb[1] * disp01[2] + pb[2] * disp01[1]) - cc * (-pc[1] * disp01[2] + pc[2] * disp01[1]);
    w4[1] = (pb[0] * disp01[2] - pb[2] * disp01[0]) - cc * (pc[0] * disp01[2] - pc[2] * disp01[0]);
    w4[2] = (-pb[0] * disp01[1] + pb[1] * disp01[0]) - cc * (-pc[0] * disp01[1] + pc[1] * disp01[0]);

    for (unsigned i = 0; i < 3; i++) {
        derivatives[0][i] = w1[i] * gamma;
        derivatives[1][i] = w4[i] * gamma;
        derivatives[2][i] = (-w1[i] - w2[i] + w3[i]) * gamma;
    }
    return true;
}

//------------------------------------------------------------
// Without derivatives.
//------------------------------------------------------------

// Calculate a squared distance.

void calc_squared_distance(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
    double rr2 = 0.0;
    std::array<double, 3> displacement;
    subtract_min_image_vectors(particle_ids, particle_positions, simulation_box_half_lengths, displacement);
    for (int i = 0; i < 3; i++) {
        rr2 += displacement[i] * displacement[i];
    }
    param_val = rr2;
}

// Calculate a distance.

void calc_distance(const std::array<int, 2> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
    calc_squared_distance(particle_ids, particle_positions, simulation_box_half_lengths, param_val);
    param_val = sqrt(param_val);
}

// Calculate the angle between three particles.

void calc_angle(const std::array<int, 3> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, double &param_val)
{   
    std::array<std::array<double, 3>, 1> dist_derivs_01;
    std::array<std::array<double, 3>, 1> dist_derivs_02;
    double rr2_01, rr2_02;
    conditionally_calc_squared_distance_and_derivatives(std::array<int, 2>({{particle_ids[0], particle_ids[1]}}), particle_positions, simulation_box_half_lengths, MAXFLOAT, rr2_01, dist_derivs_01);
    conditionally_calc_squared_distance_and_derivatives(std::array<int, 2>({{particle_ids[0], particle_ids[2]}}), particle_positions, simulation_box_half_lengths, MAXFLOAT, rr2_02, dist_derivs_02);
    
    // Calculate the cosine
    double rr_01 = sqrt(rr2_01);
    double rr_02 = sqrt(rr2_02);
    double cos_theta = dot_product(dist_derivs_01[0], dist_derivs_02[0]) / (4 * rr_01 * rr_02);
    double max = 1.0 - VERYSMALL_F;
    double min = -1.0 + VERYSMALL_F;
    if (cos_theta > max) cos_theta = max;
    else if (cos_theta < min) cos_theta = min;
    
    // Calculate the angle.
    double theta = acos(cos_theta);
    param_val = theta * DEGREES_PER_RADIAN;
}

// Calculate a dihedral angle.

void calc_dihedral(const std::array<int, 4> &particle_ids, const rvec *particle_positions, const real *simulation_box_half_lengths, double &param_val)
{
    // Find the relevant displacements for defining the angle.
    std::array<double, 3> disp20, disp01, disp23;
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[2], particle_ids[0]}}), particle_positions, simulation_box_half_lengths, disp20);
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[0], particle_ids[1]}}), particle_positions, simulation_box_half_lengths, disp01);
    subtract_min_image_vectors(std::array<int, 2>({{particle_ids[2], particle_ids[3]}}), particle_positions, simulation_box_half_lengths, disp23);
    
    // Calculate the angle, which requires many intermediates.
    double rrbc = 1.0 / sqrt(dot_product(disp01, disp01));
    std::array<double, 3> pb, pc;
    cross_product(disp20, disp01, pb);
    cross_product(disp01, disp23, pc);
    
    double pb2 = dot_product(pb, pb);
    double rpb1 = 1.0 / sqrt(pb2);
    double pc2 = dot_product(pc, pc);
    double rpc1 = 1.0 / sqrt(pc2);
    
    double pbpc = dot_product(pb, pc);
    double c = pbpc * rpb1 * rpc1;
    double s = (disp01[0] * (pc[1] * pb[2] - pc[2] * pb[1]) + disp01[1] * (pb[0] * pc[2] - pb[2] * pc[0]) + disp01[2] * (pc[0] * pb[1] - pc[1] * pb[0])) * (rpb1 * rpc1 * rrbc);
    if (s < 0.0 && s > -VERYSMALL_F) s = -VERYSMALL_F;
    if (s > 0.0 && s < VERYSMALL_F) s = VERYSMALL_F;
    param_val = atan2(s, c);
}