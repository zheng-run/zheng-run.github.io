import math

def N(a, q):
    """
    Compute N(a) = floor(a - 1) - floor((a - 1) / q).
    """
    return math.floor(a - 1) - math.floor((a - 1) / q)


def q_adic_expansion(delta_minus_1, q, max_length):
    """
    Compute the q-adic expansion of delta - 1, padded with zeros to max_length.
    Returns a list where the i-th element is delta_i.
    """
    expansion = []
    temp = delta_minus_1
    while temp > 0:
        expansion.append(temp % q)
        temp = temp // q
    # Pad with zeros to reach max_length
    while len(expansion) < max_length:
        expansion.append(0)
    return expansion


def compute_k_delta(delta_minus_1, q, h):
    """
    Compute k_delta satisfying q^{h + k_delta} <= delta - 1 < q^{h + k_delta + 1}.
    """
    k_delta = 0
    while q ** (h + k_delta) <= delta_minus_1:
        k_delta += 1
    return k_delta - 1


def find_s_delta(delta_expansion, h, k_delta, m):
    """
    Find the smallest integer s in [m - 2h - k_delta, k_delta] such that delta_{h + s} > 0.
    """
    s_range = range(m - 2 * h - k_delta, k_delta + 1)
    for s in s_range:
        if h + s < len(delta_expansion) and delta_expansion[h + s] > 0:
            return s
    return None  # If no valid s found


def compute_mu_o(delta_expansion, q, h, k_delta, s_delta):
    """
    Compute mu(delta) = min(sum_{l=0}^{h - k_delta} delta_l q^l,
    sum_{l=h - s_delta + 1}^{h + k_delta} delta_l q^{l - (h - s_delta + 1)}).
    """
    sum1 = sum(delta_expansion[l] * (q ** l) for l in range(0, h - k_delta + 1))
    sum2 = sum(delta_expansion[l] * (q ** (l - (h - s_delta + 1)))
             for l in range(h - s_delta + 1, h + k_delta + 1))
    return min(sum1, sum2)


def compute_mu_e(delta_expansion, q, h, k_delta, s_delta):
    """
    Compute \widetilde{mu}(delta) = min(sum_{l=0}^{h - k_delta - 1} delta_l q^l,
    sum_{l=h - s}^{h + k_delta} delta_l q^{l - (h - s)}).
    """
    sum1 = sum(delta_expansion[l] * (q ** l) for l in range(0, h - k_delta))
    sum2 = sum(delta_expansion[l] * (q ** (l - (h - s_delta)))
             for l in range(h - s_delta, h + k_delta + 1))
    return min(sum1, sum2)


def f_delta(delta, q, m):
    """
    Compute f(delta) (for odd m case).
    """
    h = (m - 1) // 2
    delta_minus_1 = delta - 1
    k_delta = compute_k_delta(delta_minus_1, q, h)
    delta_expansion = q_adic_expansion(delta_minus_1, q, h + k_delta + 1)

    if delta <= q ** (h + 1):
        return 0

    s_delta = find_s_delta(delta_expansion, h, k_delta, m)
    if s_delta is None:
        raise ValueError("No valid s found.")

    mu_o = compute_mu_o(delta_expansion, q, h, k_delta, s_delta)

    # Compute sum_{i=-k_delta+1}^{k_delta} sum_{t in T_i(delta)} N(t q^{2i-1} + 1)
    total_sum = 0
    for i in range(-k_delta + 1, k_delta + 1):
        lower_bound = int(q ** (k_delta - i))  # Convert to integer
        upper_bound = sum(
            delta_expansion[ell] * (q ** (ell - h - i))
            for ell in range(h + s_delta, h + k_delta + 1)
            if ell < len(delta_expansion))
        # Convert to integer
        if upper_bound == int(upper_bound):
            upper_bound -= 1
        else:
            upper_bound = int(math.floor(upper_bound))

        # Force integer conversion (to be safe)
        lower_bound = int(lower_bound)
        upper_bound = int(upper_bound)
        for t in range(lower_bound, upper_bound + 1):
            if t % q != 0:
                total_sum += N(t * (q ** (2 * i - 1)) + 1, q)

    return ((k_delta - 1) * (q ** (2 * k_delta - 3)) * ((q - 1) ** 2)
            + N(mu_o + 1, q) + total_sum)


def tilde_f_delta(delta, q, m):
    """
    Compute tilde{f}(delta) (for even m case).
    """
    h = m // 2
    delta_minus_1 = delta - 1
    k_delta = compute_k_delta(delta_minus_1, q, h)
    delta_expansion = q_adic_expansion(delta_minus_1, q, h + k_delta + 1)

    if delta <= q ** h:
        return 0
    elif delta <= q ** (h + 1):
        s_delta = find_s_delta(delta_expansion, h, k_delta, m)
        if s_delta is None:
            raise ValueError("No valid s found.")

        mu_e = compute_mu_e(delta_expansion, q, h, k_delta, s_delta)
        delta_h = delta_expansion[h]
        return 0.5 * (delta_h - 1) * delta_h + N(mu_e + 1, q)
    else:
        s_delta = find_s_delta(delta_expansion, h, k_delta, m)
        if s_delta is None:
            raise ValueError("No valid s found.")

        mu_e = compute_mu_e(delta_expansion, q, h, k_delta, s_delta)

        # Compute sum_{i=-k_delta}^{k_delta} sum_{t in T_i(delta)} N(t q^{2i} + 1)
        total_sum = 0
        for i in range(-k_delta, k_delta + 1):
            lower_bound = int(q ** (k_delta - i))  # Convert to integer
            upper_bound = sum(
                delta_expansion[ell] * (q ** (ell - h - i))
                for ell in range(h + s_delta, h + k_delta + 1)
                if ell < len(delta_expansion))  # Convert to integer
            # Convert to integer
            if upper_bound == int(upper_bound):
                upper_bound -= 1
            else:
                upper_bound = int(math.floor(upper_bound))

            # Force integer conversion (to be safe)
            lower_bound = int(lower_bound)
            upper_bound = int(upper_bound)
            for t in range(lower_bound, upper_bound + 1):
                if t % q != 0:
                    total_sum += N(t * (q ** (2 * i)) + 1, q)

        return ((k_delta - 0.5) * (q ** (2 * k_delta - 2)) * ((q - 1) ** 2)
                + 0.5 * (q ** (k_delta - 1)) * (q - 1)
                + N(mu_e + 1, q) + total_sum)


def g_delta(delta, q, m):
    """
    Compute g(delta).
    """
    h = m // 2
    delta_minus_1 = delta - 1
    delta_expansion = q_adic_expansion(delta_minus_1, q, h + 1)

    if delta <= q ** h:
        return 0
    elif delta <= q ** (h + 1):
        delta_h = delta_expansion[h]
        tau = 1 if (delta_h > 0 and
                   sum(delta_expansion[l] * (q ** l) for l in range(0, h)) >=
                   sum(delta_expansion[l] * (q ** (l - h)) for l in range(h, h + 1))) else 0
        return delta_h - 1 + tau
    else:
        k_delta = compute_k_delta(delta_minus_1, q, h)
        """
        Compute tau(delta).
        """
        tau = 1 if (delta_expansion[h] > 0 and
                   sum(delta_expansion[l] * (q ** l) for l in range(0, h)) >=
                   sum(delta_expansion[l] * (q ** (l - h)) for l in range(h, h + k_delta + 1))) else 0
        return N(sum(delta_expansion[l] * (q ** (l - h))
                for l in range(h, h + k_delta + 1)), q) + tau


def compute_dim(q, m, delta):
    """
    Compute dim(C_{(q,m,delta)}).
    """
    n = q ** m - 1
    if m % 2 == 1:  # Odd m
        f_val = f_delta(delta, q, m)
        return n - m * (N(delta, q) - f_val)
    else:  # Even m
        tilde_f_val = tilde_f_delta(delta, q, m)
        g_val = g_delta(delta, q, m)
        return n - m * (N(delta, q) - tilde_f_val) - (m // 2) * g_val


