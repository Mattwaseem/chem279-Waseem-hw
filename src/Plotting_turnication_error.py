#!/usr/bin/env python3
import matplotlib.pyplot as plt
import numpy as np

# Load the data from truncation_error_Problem2.txt
data = np.loadtxt('truncation_error_Problem2.txt', delimiter=',', skiprows=1)

# Extract h values and errors
h_values = data[:, 0]
forward_errors = data[:, 1]
central_errors = data[:, 2]

# Create the plot
plt.figure()
plt.loglog(h_values, forward_errors, 'o-', label="Forward Difference Error", marker='o', color='blue')  # Adding marker for forward difference
plt.loglog(h_values, central_errors, 'x-', label="Central Difference Error", marker='x', color='orange')  # Adding marker for central difference
plt.xlabel("Step size (h)")
plt.ylabel("Truncation Error")
plt.legend()

# Save the plot as a PNG file
plt.savefig("truncation_error_plot.png")
print("Plot saved to truncation_error_plot.png")


# Display plot
# plt.show()
