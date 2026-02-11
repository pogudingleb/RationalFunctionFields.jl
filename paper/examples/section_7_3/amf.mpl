# The code taken from the AIDA package (https://www-sop.inria.fr/members/Evelyne.Hubert/aida/) by Evelyne Hubert
# Source: https://www-sop.inria.fr/members/Evelyne.Hubert/aida/amf/amf.html

with(Groebner):

amf := proc(G::list(polynom), g::list(algebraic), P::list(polynom), z::list(name), lambda::list(name), r::symbol)
  local Z, n, Q, h, h1, R, y;


if nops(z) <> nops(g) then error "%1 should have as many components  %2", g,z
elif indets(G) minus {op(lambda)}<>{} then error "first argument (group ideal) should be polynomials in %1", lambda
elif indets(P) minus {op(z)}<>{} then error "third argument (section) should be polynomials in %1", z
elif indets(g) minus {op(lambda), op(z)} <>{} then error "second argument (group action) should be rational functions in %2", g, {op(lambda), op(z)}
end if;

n := nops(z);
Z := [ seq(Z||i, i=1..n) ];

Q := [op(G), op(map(normal,zip(`-`, Z, g) )), op(subs(zip(`=`, z, Z), P))] ;
h := convert(convert( map(p->`if`(type(p,`^`),op(1,p), p), map(denom, Q)), set), `*`);
Q := map(numer, Q);
Q := Groebner[Basis]([op(Q), h*h1-1], lexdeg([h1,op(lambda)], Z));
Q := remove(has, Q, [h1,op(lambda)]);

if P<>[] and (Q=[1] or Groebner[HilbertDimension](Q,tdeg(op(Z)))<>0) then
    error "not a cross-section to the orbits";
end if;

Q := map(`*`, map(p-> collect( p/ Groebner[LeadingCoefficient](p, tdeg(op(Z))), Z, normal, distributed), Q), -1);

R := map(coeffs, Q, Z);
R := convert(select(has, convert(map(coeffs, Q, Z), set), z), list);
if nargs=6 then
  R := [seq( R[i]=r||i, i=1..nops(R) )];
  Q := map( proc(p) local C, T;
            C := [coeffs(p, Z, T)];
            C := subs( R, C);
            convert(zip(`*`, C, [T]), `+`);
        end proc, Q);
  Q := subs( zip(`=`,Z, z), Q);
  R := map(q-> rhs(q) = lhs(q), R);
else
  Q := NULL;
end if;
map(factor,R), Q;
end proc:


