#include "LJ_PE.h"
#include <cmath>
#include <iostream>
#include <tuple>
#include <vector>

// Function to calculate the distance between two atoms
double calculate_distance(const Atom &a1, const Atom &a2)
{
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate the Lennard-Jones energy between two atoms
std::tuple<double, double, double> calculate_lj_energy_components(double distance)
{
    double ratio = sigma_au / distance;
    double ratio_6 = std::pow(ratio, 6);
    double ratio_12 = ratio_6 * ratio_6;

    double energy = 4 * epsilon_au * (ratio_12 - ratio_6);
    return {energy, ratio_6, ratio_12};
}

// Function to calculate the total Lennard-Jones energy for a cluster of atoms
double calculate_total_energy(const std::vector<Atom> &atoms)
{
    double total_energy = 0.0;
    int n = atoms.size();
    for (int i = 0; i < n; ++i)
    {
        for (int j = i + 1; j < n; ++j)
        {
            double distance = calculate_distance(atoms[i], atoms[j]);

            std::cout << "Distance between atom " << i << " and atom " << j << ": " << distance << " Å\n";
            if (distance < 0.01)
            {
                std::cout << "Warning: Very small distance between atom " << i << " and atom " << j << ": " << distance << " Å\n";
                // There's no need to add the (distance < 0.01) statement again here; it's already inside the if block.
            }

            // Reuse the energy function and extract only the energy
            auto [lj_energy, sigma_r_6, sigma_r_12] = calculate_lj_energy_components(distance);

            std::cout << "Energy contribution from atom pair " << i << "-" << j << ": " << lj_energy << " kcal/mol\n";

            total_energy += lj_energy;
            std::cout << "Total energy so far: " << total_energy << " kcal/mol\n";
        }
    }
    return total_energy;
}

// Function to calculate the force between two atoms
std::vector<double> calculate_lj_force(const Atom &a1, const Atom &a2)
{
    std::vector<double> force(3, 0.0); // Force vector [Fx, Fy, Fz]

    // Reuse the distance calculation from Problem 1
    double r = calculate_distance(a1, a2); // Distance between atoms

    // Reuse the energy function to get ratios
    auto [energy, sigma_r_6, sigma_r_12] = calculate_lj_energy_components(r);

    // Calculate the magnitude of the force
    double force_magnitude = 24 * epsilon_au * (2 * sigma_r_12 / (r * r) - sigma_r_6 / (r * r));

    // Now apply this force along the direction of r_ij (normalized vector)
    force[0] = force_magnitude * ((a1.x - a2.x) / r); // Fx
    force[1] = force_magnitude * ((a1.y - a2.y) / r); // Fy
    force[2] = force_magnitude * ((a1.z - a2.z) / r); // Fz

    return force;
}

// calculating forward diff appropximations
double forward_difference_force(const std::vector<Atom> &atoms, int atom_index, double h)
{
    // Store the original position
    Atom original_atom = atoms[atom_index];

    // Modify position by +h in the x-direction
    std::vector<Atom> atoms_plus_h = atoms;
    atoms_plus_h[atom_index].x += h;

    // Calculate energy at R + h
    double energy_plus_h = calculate_total_energy(atoms_plus_h);

    // Calculate energy at the original position
    double energy_original = calculate_total_energy(atoms);

    // Forward difference approximation for the force
    return -(energy_plus_h - energy_original) / h;
}

//  calculating central diff approx
double central_difference_force(const std::vector<Atom> &atoms, int atom_index, double h)
{
    // Store original position
    Atom original_atom = atoms[atom_index];

    // Modify position by +h in the x-direction
    std::vector<Atom> atoms_plus_h = atoms;
    atoms_plus_h[atom_index].x += h;

    // Modify position by -h in the x-direction
    std::vector<Atom> atoms_minus_h = atoms;
    atoms_minus_h[atom_index].x -= h;

    // Calculate energy at R + h
    double energy_plus_h = calculate_total_energy(atoms_plus_h);

    // Calculate energy at R - h
    double energy_minus_h = calculate_total_energy(atoms_minus_h);

    // Central difference approximation for the force
    return -(energy_plus_h - energy_minus_h) / (2 * h);
}
