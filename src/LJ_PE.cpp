#include "LJ_PE.h"
#include <cmath>
#include <iostream>

// Function to calculate the distance between two atoms
double calculate_distance(const Atom &a1, const Atom &a2)
{
    double dx = a1.x - a2.x;
    double dy = a1.y - a2.y;
    double dz = a1.z - a2.z;
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

// Function to calculate the Lennard-Jones energy between two atoms
double calculate_lj_energy(double distance)
{
    double ratio = sigma_au / distance;
    double ratio_6 = std::pow(ratio, 6);
    double ratio_12 = ratio_6 * ratio_6;
    return 4 * epsilon_au * (ratio_12 - ratio_6);
}

// Function to calculate the total Lennard-Jones energy for a cluster of atoms
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
            }

            double lj_energy = calculate_lj_energy(distance);
            std::cout << "Energy contribution from atom pair " << i << "-" << j << ": " << lj_energy << " kcal/mol\n";

            total_energy += lj_energy; // Add only once
        }
    }
    return total_energy;
}
