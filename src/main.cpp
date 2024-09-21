#include <iostream>
#include <filesystem>
#include "file_io.h"
#include "LJ_PE.h"

int main()
{
    try
    {
        std::string directory_path = "sample_input/Energy";
        double h = 0.01; // Step size for finite difference

        // Loop through all files in the specified directory
        for (const auto &entry : std::filesystem::directory_iterator(directory_path))
        {
            std::string filename = entry.path().string();
            std::vector<Atom> atoms = read_atoms_from_file(filename);
            std::cout << "Processing file: " << filename << "\n";
            std::cout << "Atoms in the cluster:\n";

            for (const auto &atom : atoms)
            {
                std::cout << atom.atomic_number << " " << atom.x << " " << atom.y << " " << atom.z << "\n";
            }

            // Calculate total Lennard-Jones energy
            double total_energy = calculate_total_energy(atoms);
            std::cout << "Total Lennard-Jones energy: " << total_energy << " kcal/mol\n";

            // Compare forces for each atom (analytical vs forward difference vs central difference)
            for (int i = 0; i < atoms.size(); ++i)
            {
                std::cout << "\nComparing forces for atom " << i << ":\n";

                // Analytical force between atom 0 and atom i
                std::vector<double> analytical_force = calculate_lj_force(atoms[0], atoms[i]);

                // Forward difference force for atom i
                double forward_force = forward_difference_force(atoms, i, h);

                // Central difference force for atom i
                double central_force = central_difference_force(atoms, i, h);

                // Output results
                std::cout << "  Analytical force: Fx = " << analytical_force[0]
                          << ", Fy = " << analytical_force[1]
                          << ", Fz = " << analytical_force[2] << "\n";
                std::cout << "  Forward difference force (Fx): " << forward_force << "\n";
                std::cout << "  Central difference force (Fx): " << central_force << "\n";
            }
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}
