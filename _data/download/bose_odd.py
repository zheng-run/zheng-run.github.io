import math

# m is odd

def q_adic_expansion(delta, q):
    """
    Compute the q-adic expansion of delta.
    Returns a list where the i-th element is delta_i.
    """
    expansion = []
    while delta > 0:
        expansion.append(delta % q)
        delta = delta // q
    return expansion


def compute_k_delta(delta, q, h):
    """
    Compute k_delta satisfying q^{h + k_delta} <= delta < q^{h + k_delta + 1}.
    """
    k_delta = 0
    while q ** (h + k_delta) <= delta:
        k_delta += 1
    return k_delta - 1


def find_s(delta_expansion, h, k_delta):
    """
    Find the smallest integer s in [-k_delta + 1, k_delta] such that delta_{h + s} > 0.
    """
    for s in range(-k_delta + 1, k_delta + 1):
        if delta_expansion[h + s] > 0:
            return s
    return None


def compute_hat_delta(delta_expansion, q, h, s, k_delta):
    """
    Compute hat_delta = sum_{l=h+s}^{h+k_delta} delta_l q^l + sum_{l=h-s+1}^{h+k_delta} delta_l q^{l - (h - s + 1)}.
    """
    hat_delta = 0
    # First part: sum_{l=h+s}^{h+k_delta} delta_l q^l
    for l in range(h + s, h + k_delta + 1):
        hat_delta += delta_expansion[l] * (q ** l)
    # Second part: sum_{l=h-s+1}^{h+k_delta} delta_l q^{l - (h - s + 1)}
    for l in range(h - s + 1, h + k_delta + 1):
        hat_delta += delta_expansion[l] * (q ** (l - (h - s + 1)))
    return hat_delta


def compute_d_B(q, m, delta):
    """
    Compute the value of d_B(delta).
    """
    h = (m - 1) // 2  # h = (m - 1) / 2
    k_delta = compute_k_delta(delta, q, h)
    delta_expansion = q_adic_expansion(delta, q)
    s = find_s(delta_expansion, h, k_delta)

    if s is None:
        raise ValueError("No valid s found.")

    hat_delta = compute_hat_delta(delta_expansion, q, h, s, k_delta)

    # Compute sum_{l=0}^{h - k_delta} delta_l q^l
    sum_left = sum(delta_expansion[l] * (q ** l) for l in range(0, h - k_delta + 1))
    # Compute sum_{l=h - s + 1}^{h + k_delta} delta_l q^{l - (h - s + 1)}
    sum_right = sum(delta_expansion[l] * (q ** (l - (h - s + 1))) for l in range(h - s + 1, h + k_delta + 1))

    if delta % q != 0:  # q does not divide delta
        if sum_left > sum_right:
            return delta
        elif delta_expansion[h - s + 1] != q - 1:
            return hat_delta + 1
        else:
            return hat_delta + 2
    else:  # q divides delta
        if sum_left >= sum_right:
            return delta + 1
        elif delta_expansion[h - s + 1] != q - 1:
            return hat_delta + 1
        else:
            return hat_delta + 2


# Coset leader method verification
def q_cyclotomic_cosets(q, m):
    """
    Compute all q-cyclotomic cosets modulo q^m - 1.
    Returns a dictionary where keys are coset leaders and values are the cosets.
    """
    n = q ** m - 1
    visited = set()
    cosets = {}

    for a in range(n + 1):
        if a not in visited:
            # Generate the coset for a
            coset = []
            current = a
            while current not in coset:
                coset.append(current)
                visited.add(current)
                current = (current * q) % n
            # The smallest element in the coset is the coset leader
            coset_leader = min(coset)
            cosets[coset_leader] = coset
    return cosets


def find_next_coset_leader(cosets, delta):
    """
    Find the smallest coset leader greater than delta - 1.
    """
    # Get all coset leaders
    leaders = sorted(cosets.keys())

    # Find the smallest leader > delta - 1
    for leader in leaders:
        if leader > delta - 1:
            return leader
    return None  # If no such leader exists


def coset_compute_dimension_and_dB(q, m, delta):
    """
    Compute the dimension of the code C_{(q,m,delta,b)} and d_B(delta).
    d_B(delta) is the smallest coset leader greater than delta - 1.
    """
    n = q ** m - 1
    cosets = q_cyclotomic_cosets(q, m)

    # Filter coset leaders in [1, delta - 1]
    relevant_coset_leaders = [a for a in cosets.keys() if 1 <= a <= delta - 1]

    # Sum the sizes of the relevant cosets
    total_size = sum(len(cosets[a]) for a in relevant_coset_leaders)

    # Compute the dimension
    dimension = n - total_size

    # Find d_B(delta)
    d_B = find_next_coset_leader(cosets, delta)
    if d_B is None:
        raise ValueError(f"No coset leader greater than {delta - 1} exists.")

    return dimension, d_B

