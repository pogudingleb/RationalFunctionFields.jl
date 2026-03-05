read "amf.mpl";

with(LinearAlgebra):

lambda := [seq(seq(v[i, j], i=1..3), j=1..3)];
z := [a1, a2, a3, b1, b2, b3];

G := [
    seq(seq(v[i, j] * (1 - v[i, j]), i=1..3), j=1..3),
    seq(v[i, 1] + v[i, 2] + v[i, 3] - 1, i=1..3),
    seq(v[1, i] + v[2, i] + v[3, i] - 1, i=1..3)
];

V := Matrix([seq([seq(v[i, j], i=1..3)], j=1..3)]);

d := Determinant(V);

gg := [
  op(convert(V . Matrix([[a1], [a2], [a3]]), list)),
  op(convert(V . Matrix([[d * b1], [d * b2], [d * b3]]), list))
];


invar := amf( G, gg, [], z, lambda);

with(FileTools[Text]):

fname := "example_7_9.gens";
for i from 1 to nops(invar) do
    WriteString(fname, convert(invar[i], string));
    WriteString(fname, ",\n");
end do;
Close(fname);
