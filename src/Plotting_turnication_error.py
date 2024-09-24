import matplotlib.pyplot as plt
import numpy as np

# Load the data from the CSV file generated from main.cpp file 
data = np.loadtxt('truncation_error_Problem2.txt', delimiter=',', skiprows=1)

h = data[:, 0]
forward_error = data[:, 1]
central_error = data[:, 2]

# Plotting the log-log plot
plt.loglog(h, forward_error, label='Forward Difference Error', marker='o')
plt.loglog(h, central_error, label='Central Difference Error', marker='s')

plt.xlabel('Step size (h)')
plt.ylabel('Error')
plt.title('Truncation Error vs Step Size')
plt.legend()
plt.grid(True)

# Show the plot
plt.show()
