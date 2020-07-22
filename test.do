clear 
discard 

adopath ++ "src/"

sysuse auto 

egen dummy = cut(headroom), group(3)

stepdown (reg price mpg, adjust(mpg) ) /// 
    (reg price weight, vce(cluster foreign) adjust(weight) ), reps(20)

mat list r(p)
mat list r(adj_p)

di "coeflist | $sdCoefList"