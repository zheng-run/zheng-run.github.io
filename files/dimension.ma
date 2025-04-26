q := 3;  
m := 4;  
n := q^m-1; // Example: length of the code


F := GF(q);

// Calculate delta_max
delta_max := q^(Floor((2*m - 1)/3) + 1);

// Loop through delta values
for delta in [2..delta_max] do
    BCHCode := BCHCode(F, n, delta);
    dimension := Dimension(BCHCode);
    print  dimension;
end for;