# stepdown
 flexibly implements the Westfall-Young step-down resampling method 

existing packages are powerful but have difficult to use syntax. This implementation has very simple user-facing syntax. 

# syntax

inspired by `twoway`, this command allows for a list of estimation commands with the extra option `adjust`. 

We call `stepdown` by typing sets of estimatation commands, enclosed in parenthesis. We then flag the coefficients with p-values we'd like to adjust with the `adjust` options. There are no restrictions requiring the estimation commands to be the same, or for the arguments to `adjust` to be the same. The FWER will be controlled over all coefficients flagged in a particular call to stepdown. In this example, we will control the FWER over estimates of coeffients on `x2` and `x3` in the first model and coefficients on `x2` and `x4` in the second model. 

``` 
stepdown (reg y1 x1 x2 x3, adjust(x2 x3)) (areg y2 x2 x4 x5,  absorb(z) adjust(x2 x4))
```

Except for the `adjust` option, estimation commands will be submitted exactly as specified here. Any suboptions supported by the underlying estimation command are allowed, like clustering or high dimensional fixed effects with `reghdfe`. Additionally, `adjust` supports factor variable syntax. So, the coefficients for particular factor estimates can also be added to the list of coefficients to adjust. 

``` 
stepdown (reg y1 i.z, adjust(1.z 2.z 3.z)) (areg y2 x2 x4 x5,  absorb(z) adjust(x2 x4))
```

# meta-options 

Since `stepdown` uses a bootstrap procedure, `stepdown` allows bootstrap options to be appended to the end of the function call. 

In particular, `stepdown` recognizes `reps` to control the number of bootstrap repitions (1000 is recommended, but the default is only 10). `strata` `cluster` `idcluster` and `weight` will all be passed through to `gsample` (see `help gsample`). Additionally, `bs_sample` takes the bootstrap sample size (defaults to a bootstrap draw of equivalent size). By default, `stepdown` runs quietly, but the `echo` option will cause `stepdown` to announce each bootstrap iteration. 

``` 
stepdown (reg y1 x1 x2 x3, adjust(x2 x3)) (areg y2 x2 x4 x5,  absorb(z) adjust(x2 x4)), reps(1000) cluster(classroom) strata(state)
```

# return

`stepdown` returns two matrices to `r()`. It returns `r(p)` which a list of the original, unadjusted p-values, and `r(adj_p)` which contains the Westfall-Young adjusted p-values. 
