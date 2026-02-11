read "amf.mpl";

# Code for the example is taken from https://www-sop.inria.fr/members/Evelyne.Hubert/aida/amf/amf.html

lambda := [seq(lambda||i, i=1..4)];
z := [seq(z||i, i=1..7)];

G := [lambda1 * lambda4-lambda3 * lambda2 - 1];

g := [
    lambda1 * z1 + lambda2 * z2, 
    lambda3 * z1 + lambda4 * z2,
    lambda1 * z3 + lambda2 * z4, 
    lambda3 * z3 + lambda4 * z4,
    z7 * lambda2^2 + 2 * z6 * lambda2 * lambda1 + z5 * lambda1^2,
    z5 * lambda3 * lambda1 + (1 + 2 * lambda3 * lambda2) * z6 + z7 * lambda4 * lambda2,
    z5 * lambda3^2 + z7 * lambda4^2 + 2 * z6 * lambda4 * lambda3
];

inv := amf(G, g, [], z, lambda);
lprint(inv);
