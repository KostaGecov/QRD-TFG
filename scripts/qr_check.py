import numpy
import os

numpy.set_printoptions(precision=4, linewidth=300)

A = numpy.loadtxt("data_files/data_in_16.dat", dtype=float)

Q, R = numpy.linalg.qr(A)

Q = numpy.array(numpy.around(Q, decimals=26))
R = numpy.array(numpy.around(R, decimals=26))
A_res = numpy.dot(Q, R)

print("\nA:\n", A)
print("\nQ:\n", Q)
print("\nR:\n", R)
print("\nA_res:\n", A_res)

# Define the new directory
new_directory = "data_files"

# Create the new directory if it doesn't exist
if not os.path.exists(new_directory):
    os.makedirs(new_directory)

# Specify the full path to save the file under the data_files directory
file_path = os.path.join(new_directory, "data_out_gold_16.dat")

with open(file_path, "w") as dataOut:
    for row in R:
        dataOut.write(" ".join([str(a) for a in row]) + "\n")
