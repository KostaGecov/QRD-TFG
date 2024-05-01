import numpy as np
import os

# Define the dimensions of the matrix
rows = 16
cols = 16

# Define the interval for random numbers
lower_bound = -100
upper_bound = 100

# Generate a dense random matrix
dense_random_matrix = np.random.uniform(lower_bound, upper_bound, size=(rows, cols))

# Define the new directory
new_directory = "data_files"

# Create the new directory if it doesn't exist
if not os.path.exists(new_directory):
    os.makedirs(new_directory)

# Specify the full path to save the file under the data_files directory
file_path = os.path.join(new_directory, "data_in_16.dat")

# Save the matrix to a .dat file
np.savetxt(file_path, dense_random_matrix)

print(f"Matrix saved to {file_path}")
