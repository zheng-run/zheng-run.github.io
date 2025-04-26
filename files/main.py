import math
import dim
import bose_odd as bo
import bose_even as be

# Given q and m, compute the dimension of C(q, m, delta) for delta  [2, q^(floor((2m-1)/3) + 1)]
# and compute the Bose distance of C(q, m, delta) for delta [2, q^(floor((2m-1)/3) + 1)).

q = 3
m=7
h = math.floor(m / 2)
n = q**m - 1
upper_bound = int(q ** (math.floor((2 * m - 1) / 3) + 1))

for delta in range(2, upper_bound + 1):
    # Compute the dimension of the code C(q, m, delta)
    k = dim.compute_dim(q, m, delta)
    k = int(k)
    # Compute the Bose distance (dB) based on the value of delta
    if delta <= q**(m - h) - 1:
        # Case 1: delta is small (<= q^{m-h} - 1)
        if delta % q == 0:
            dB = delta + 1  # Adjust if delta is divisible by q
        else:
            dB = delta  # Otherwise, dB = delta
    elif q**(m - h) <= delta < upper_bound:
        # Case 2: delta is in the middle range
        if m % 2 == 0:
            dB = be.compute_dB(q, m, delta)  # Use even-m Bose distance formula
        else:
            dB = bo.compute_d_B(q, m, delta)  # Use odd-m Bose distance formula
    else:
        # Case 3: delta = q^(floor((2m-1)/3) + 1)
        dB = None
    # Print results
    print("Bose distance of C(",q,m,delta,")=", dB,";Dimension of C(",q,m,delta,")=", k)

