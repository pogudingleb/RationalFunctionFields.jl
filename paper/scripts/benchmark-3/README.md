# Simplification in Julia and Maple

Script for running Julia implementation of our algorithm and the Maple implementation of the one from (https://arxiv.org/abs/2004.07774) on our benchmark suite and collecting various statistics.

## In Julia

1. Generate scripts:

```
julia populate.jl
```

2. Run some of the scripts (say, for only some models):

- Run our simplification algorithm for Goodwin and SLIQR:
```
python ../run.py --pattern='benchmark-3 & simplify & Goodwin.jl | SLIQR.jl' --timeout=3600
```

- Collect statistics about the original generating sets for Goodwin and SLIQR:
```
python src/run.py --pattern='benchmark-3 & input_stats & Goodwin.jl | SLIQR.jl" --timeout=3600
```

3. Collect results:

```
julia collect.jl
```

## In Maple

0. Get `AllIdentifiableFunctions`:

```
git clone https://github.com/pogudingleb/AllIdentifiableFunctions
```

1. Same as in Julia.

2. First, generate the input generating sets for Maple:

```
python ../run.py --pattern='benchmark-3 & maple_simplify & generate_input & Goodwin | SLIQR' --timeout=3600
```

Then, run Maple code with:

```
python ../run.py --pattern='benchmark-3 & maple_simplify & Goodwin.mpl | SLIQR.mpl' --timeout=3600
```

3. Same as in Julia.
