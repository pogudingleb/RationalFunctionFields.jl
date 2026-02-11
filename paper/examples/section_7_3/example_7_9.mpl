read "amf.mpl";

with(LinearAlgebra):

M := Matrix([[a, b, c], [0, a_inv, 0], [0, 0, a^2]]):
z := [seq(z||i, i=1..9)];
Mz := Matrix([[z[1], z[2], z[3]], [z[4], z[5], z[6]], [z[7], z[8], z[9]]]):

G := [a * a_inv - 1];

A := M . Mz;
g := [ seq(seq(A[i, j], j = 1..3), i=1..3)];

inv := amf(G, g, [], z, [a, a_inv, b, c]);
lprint(inv);
