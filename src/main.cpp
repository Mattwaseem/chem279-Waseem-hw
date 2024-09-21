#include <iostream>
#include <filesystem>
#include "file_io.h"
#include "LJ_PE.h"

int main()
{
    try
    {
        std::string directory_path = "sample_input/Energy";

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

            double total_energy = calculate_total_energy(atoms);
            std::cout << "Total Lennard-Jones energy: " << total_energy << " kcal/mol\n";
        }
    }
    catch (const std::exception &e)
    {
        std::cerr << e.what() << '\n';
        return 1;
    }

    return 0;
}
