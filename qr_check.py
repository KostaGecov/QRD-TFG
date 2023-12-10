import numpy

numpy.set_printoptions(precision=4, linewidth=300)

# RandomMatrix = numpy.random.randint(-1, 1, (256, 256))
A = numpy.random.uniform(-1, 1, size=(12, 12))
# print('\nrandomMatrix:\n', A)

# numpy.savetxt('data_in2.dat', randomMatrix)

# A = numpy.loadtxt("data_in2.dat", dtype=float)


# A = numpy.array(numpy.around(A, decimals=15))
Q, R = numpy.linalg.qr(A)

Q = numpy.array(numpy.around(Q, decimals=15))
R = numpy.array(numpy.around(R, decimals=15))
A_res = R*Q

print("\nA:\n", A)
print("\nQ:\n", Q)
print("\nR:\n", R)
print("\nA_res:\n", A_res)

with open("data_out.dat", "w") as dataOut:
    for row in R:
        dataOut.write(" ".join([str(a) for a in row]) + "\n")

# A_Reconstructed = numpy.around((numpy.dot(Q, R_Vitis)), decimals = 15)
# print('\nA_Reconstructed:\n', A_Reconstructed)
