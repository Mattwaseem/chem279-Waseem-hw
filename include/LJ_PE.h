#ifndef LJ_PE_H
#define LJ_PE_H

#include <vector>
#include "file_io.h" // Include for the Atom structure

// Constants for gold (Au) atoms
const double epsilon_au = 5.29; // kcal/mol
const double sigma_au = 2.951;  // Angstrom

// Function declarations for Lennard-Jones energy and distance calculation
double calculate_distance(const Atom &a1, const Atom &a2);
double calculate_lj_energy(double distance);
double calculate_total_energy(const std::vector<Atom> &atoms);

#endif // LJ_PE_H
