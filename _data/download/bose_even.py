def compute_k_delta(q, h, delta):
    """
    Compute k_delta such that q^(h + k_delta) <= delta < q^(h + k_delta + 1).

    Args:
        q (int): Base number.
        h (int): h = m / 2 (half of m).
        delta (int): Integer in [q^h, q^{floor((2m-1)/3) + 1}).

    Returns:
        int: Value of k_delta.
    """
    k_delta = 0
    while q ** (h + k_delta) <= delta:
        k_delta += 1
    return k_delta - 1


def q_adic_expansion(q, delta, length):
    """
    Compute the q-adic expansion of delta.

    Args:
        q (int): Base number.
        delta (int): Number to expand.
        length (int): Length of the expansion.

    Returns:
        list: q-adic expansion of delta (LSB first).
    """
    expansion = []
    for _ in range(length):
        expansion.append(delta % q)
        delta = delta // q
    return expansion


def compute_s(h, k_delta, delta_expansion):
    """
    Find the smallest s ∈ [-k_delta, k_delta] where delta_{h+s} > 0.

    Args:
        h (int): h = m / 2.
        k_delta (int): Value of k_delta.
        delta_expansion (list): q-adic expansion of delta.

    Returns:
        int: Value of s, or None if not found.
    """
    for s in range(-k_delta, k_delta + 1):
        if delta_expansion[h + s] > 0:
            return s
    return None


def compute_hat_delta(q, h, s, k_delta, delta_expansion):
    """
    Compute hat_delta value.

    Args:
        q (int): Base number.
        h (int): h = m / 2.
        s (int): Value of s.
        k_delta (int): Value of k_delta.
        delta_expansion (list): q-adic expansion of delta.

    Returns:
        int: Value of hat_delta.
    """
    hat_delta = 0
    for ell in range(h + s, h + k_delta + 1):
        hat_delta += delta_expansion[ell] * (q ** ell)
    for ell in range(h - s, h + k_delta + 1):
        hat_delta += delta_expansion[ell] * (q ** (ell - (h - s)))
    return hat_delta


def compute_dB(q, m, delta):
    """
    Compute the Bose distance d_B(delta).

    Args:
        q (int): Base number.
        m (int): Even integer >= 4.
        delta (int): Integer in [q^h, q^{floor((2m-1)/3) + 1}).

    Returns:
        int: Value of d_B(delta).

    Raises:
        ValueError: If no valid s is found.
    """
    h = m // 2

    # Compute k_delta
    k_delta = compute_k_delta(q, h, delta)

    # Compute q-adic expansion of delta
    delta_expansion = q_adic_expansion(q, delta, h + k_delta + 1)

    # Find s
    s = compute_s(h, k_delta, delta_expansion)

    if s is None:
        raise ValueError("No valid s found.")

    # Compute hat_delta
    hat_delta = compute_hat_delta(q, h, s, k_delta, delta_expansion)

    # Compute sum1 and sum2
    sum1 = sum(delta_expansion[ell] * (q ** ell) for ell in range(0, h - k_delta))
    sum2 = sum(delta_expansion[ell] * (q ** (ell - (h - s))) for ell in range(h - s, h + k_delta + 1))

    # Calculate d_B(delta) based on conditions
    if s != 0:
        if delta % q != 0:  # q doesn't divide delta
            if sum1 > sum2:
                dB = delta
            elif delta_expansion[h - s] != q - 1:
                dB = hat_delta + 1
            else:
                dB = hat_delta + 2
        else:  # q divides delta
            if sum1 >= sum2:
                dB = delta + 1
            elif delta_expansion[h - s] != q - 1:
                dB = hat_delta + 1
            else:
                dB = hat_delta + 2
    else:  # s == 0
        if sum1 >= sum2:
            dB = delta + 1 if delta % q == 0 else delta
        else:
            dB = hat_delta

    return dB