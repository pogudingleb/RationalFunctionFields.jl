read "amf.mpl";

with(LinearAlgebra):

lambda := [lambda1, lambda2, lambda3, lambda4, c1];
z := [a10, a01, a20, a11, a02, b10, b01, b20, b11, b02];

G := [
    lambda1 * lambda4 - 1, c1 * (c1 - 1)
];

V := Matrix([
    [lambda1, lambda2, 0,         0,                     0                ],
    [0,       lambda4, 0,         0,                     0                ],
    [0,       0,       lambda1^2, 2 * lambda1 * lambda2, lambda2^2        ],
    [0,       0,       0,         lambda4 * lambda1,     lambda4 * lambda2],
    [0,       0,       0,         0,                     lambda4^2        ]
]);

A := Matrix([[a10], [a01], [a20], [a11], [a02]]);
B := Matrix([[b10], [b01], [b20], [b11], [b02]]);

gg := [
  op(convert(V . (c1 * A + (1 - c1) * B), list)),
  op(convert(V . ((1 - c1) * A + c1 * B), list))
];

invar := amf( G, gg, [a01 - 1, b10], z, lambda);


with(FileTools[Text]):

fname := "example_7_10.gens";
for i from 1 to nops(invar) do
    WriteString(fname, convert(invar[i], string));
    WriteString(fname, ",\n");
end do;
Close(fname);
