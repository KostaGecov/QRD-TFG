import numpy

numpy.set_printoptions(precision=4, linewidth=300)

A = numpy.loadtxt("data_in.dat", dtype=float)

Q, R = numpy.linalg.qr(A)

Q = numpy.array(numpy.around(Q, decimals=26))
R = numpy.array(numpy.around(R, decimals=26))
A_res = numpy.dot(Q, R)

print("\nA:\n", A)
print("\nQ:\n", Q)
print("\nR:\n", R)
print("\nA_res:\n", A_res)

with open("data_out_gold.dat", "w") as dataOut:
    for row in R:
        dataOut.write(" ".join([str(a) for a in row]) + "\n")
