read "amf.mpl";

GetM := proc(i, j)
    convert(cat("m_", i, "_", j), symbol):
end proc;

group_eqs := [a^2 + b^2 - 1];
action := [];
ms := [];
DEG := 3;

for i from 0 to DEG do
    for j from 0 to DEG - i do
        if i + j > 0 then
            exprs := expand((a * x + b * y)^i * (-b * x + a * y)^j);
            res := add(seq(seq(GetM(s, t) * coeff(coeff(exprs, x, s), y, t), t=0..DEG), s=0..DEG));
            ms := [op(ms), GetM(i, j)];
            action := [op(action), res];
        end if:
    end do:
end do:

inv := amf(group_eqs, action, [], ms, [a, b]);
lprint(inv);
