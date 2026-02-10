## Benchmarks

Benchmarking implementations for computing Groebner bases in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

#### How to run

1. To generate directories with benchmark scripts, run:

```julia
julia populate.jl "Bilirubin,EAIHRD,SLIQR,Lincomp2,Pharm,MAPK_5out,Fujita,Goodwin,i4_boku,param1,param2,simson3"
```

2. The generated scripts can be run individually. Alternatively, use the following command to run many at a time:

```
python ../run.py --pattern='benchmark-1, Fujita | Goodwin'
```

(check `python ../run.py --help` for options, such as specifying the path for Maple and Singular)

#### Benchmarked methods

1. f4-direct. Groebner.jl, directly computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$.

2. f4-flat. Groebner.jl, computing in $\mathbb{Q}[x_1,\ldots,x_m, y_1,\ldots,y_n]$ using a block ordering with $x_1,\ldots,x_m < y_1,\ldots,y_n$.

3. ffmodstd. ffmodstd package in Singular, uses interpolation.

4. paramgb. ParamPunPam.jl, our algorithm, computing in $\mathbb{Q}(x_1,\ldots,x_m)[y_1,\ldots,y_n]$ using sparse interpolation.

5. slimgb. slimgb package in Singular.
