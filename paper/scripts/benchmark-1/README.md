## Benchmarks

31 December, 2025. Alexander Demin.

Benchmarking the computation of Groebner bases in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

#### How to run

1. Do `julia populate.jl` to generate directories with benchmark scripts (see below). The generated scripts can be run individually. 

#### Benchmark scripts

> f4-direct

Method: Groebner.jl, directly computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

> f4-flat

Method: Groebner.jl, computing in $\mathbb{Q}[x_1,\ldots,x_m, y_1,\ldots,y_n]$ using a block ordering with $x_1,\ldots,x_m < y_1,\ldots,y_n$.

> ffmodstd

Method: ffmodstd package in Singular, uses interpolation.

> paramgb

Method: ParamPunPam.jl, computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$ using sparse interpolation.

> slimgb

Method: slimgb package in Singular.
