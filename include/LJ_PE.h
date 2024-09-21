#ifndef LJ_PE_H
#define LJ_PE_H

#include <vector>
#include "file_io.h" // Include for the Atom structure

// Constants for gold (Au) atoms
const double epsilon_au = 5.29; // kcal/mol
const double sigma_au = 2.951;  // Angstrom

// Function declarations for Lennard-Jones energy and distance calculation
double calculate_distance(const Atom &a1, const Atom &a2);
std::tuple<double, double, double> calculate_lj_energy_components(double distance);
double calculate_total_energy(const std::vector<Atom> &atoms);

// Function to declareand calculate force
std::vector<double> calculate_lj_force(const Atom &a1, const Atom &a2);

// forward diff and central diff eqxn fxns 
double forward_difference_force(const std::vector<Atom> &atoms, int atom_index, double h);
double central_difference_force(const std::vector<Atom> &atoms, int atom_index, double h);

#endif // LJ_PE_H
