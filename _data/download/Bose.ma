// Define a function to calculate the q-cyclotomic coset leader
CyclotomicCosetLeader := function(q, m, delta)
    n := q^m - 1; // Calculate n
    leaders := []; // Initialize an empty list to store coset leaders

    // Iterate over all possible a
    for a in [0..n-1] do
        // Calculate coset C_a
        coset := [];
        k := 0;
        repeat
            element := (a * q^k) mod n;
            if not element in coset then
                Append(~coset, element);
                k := k + 1;
            else
                break; // Stop if the element is already in the coset
            end if;
        until false;

        // Find the coset leader
        coset_leader := Minimum(coset);
        Append(~leaders, coset_leader);
    end for;

    // Find the smallest coset leader greater than or equal to delta
    min_leader := Infinity();
    for leader in leaders do
        if leader ge delta and leader lt min_leader then
            min_leader := leader;
        end if;
    end for;

    // Return the result
    if min_leader eq Infinity() then
        return "No coset leader found greater than or equal to delta";
    else
        return min_leader;
    end if;
end function;

// Example usage
q := 3; // Choose q
m := 5; // Choose m

// Calculate the upper limit for delta
delta_max := q^(Floor((2*m - 1)/3)+1);

// Iterate over delta values from 2 to delta_max
for delta in [2..delta_max-1] do
    result := CyclotomicCosetLeader(q, m, delta);
    print "the Bose distance of C(",q,m, delta, ")is:", result;
end for;