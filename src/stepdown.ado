
program stepdown, rclass
set seed 3424421

syntax anything(id="model list" name=commands everything), ///
    [Reps(integer 3)] /// /* set number of bootstrap draws  */
    [BS_sample(integer -1)] /// /* set bootstrap sample size  */
    [echo] /// option to control if it should display every bootstrap iteration 
    [STRata(passthru) CLuster(passthru) IDcluster(passthru) Weight(passthru)] /* passthru args for bootstrap */


/* default to bootstrapping a sample of the same size as the original data*/
if `bs_sample' == -1 local bs_sample

tempname pmat adj_mat model_names sdCoefIds bs_pvals
global sdPVal `pmat'
global sdPAdjusted `adj_mat'
global sdCoefIds `sdCoefIds'
global sdBootPVals `bs_pvals'

global sdCoefList 
global sdModelIndex 1


local input_command_list `commands'

/* get the initial list of p-values */
while "`commands'" != ""{
    
    gettoken cmd commands : commands, parse(" ||") bind match(parenthesis)
    
    /* di as text "Executing command: " as input `"`cmd'"' */
    
    qui: store_pvalues `cmd' mat($sdPVal)

    /* keep track of the number of models */
    global sdModelIndex = ${sdModelIndex}+1
}



/* make sure we can label the matrix with the right names after we sort them */
local n_rows : list sizeof global(sdCoefList)
matrix $sdCoefIds = 1 
forvalues i = 2/`n_rows'{
    matrix $sdCoefIds = $sdCoefIds \ `i'
}


/* sort the p-values to be ascending */
mata p = st_matrix("$sdPVal")
mata ids = st_matrix("$sdCoefIds")
mata ordered = sort((p,ids), 1)
mata st_matrix("$sdPVal",  ordered[,1])
mata st_matrix("$sdCoefIds",  ordered[,2])
 
 /* relabel now that we've sorted  */
local sorted_coef_list 
forvalues i = 1/`n_rows'{
    local id = $sdCoefIds[`i',1]
    local name :word `id' of $sdCoefList
    local sorted_coef_list `sorted_coef_list' `=strtrim("`name'")'
}  

matrix rownames $sdPVal  = `sorted_coef_list'

******************** BEGIN BOOTSTRAP *****************


/* now we bootstrap and re-estimate the coefficient vector */
forvalues bs_index = 1/`reps'{
    preserve 

    if !missing("`echo'") di as text "bootstrap iteration `bs_index'"
    // we support all the passthrough args for bsample 
    if !missing("`strata'`cluster'`idcluster'`weight'") local comma ","
    bsample `bs_sample' `comma' `strata' `cluster' `idcluster' `weight' 

    /* get the original list of commands to loop through */
    local commands `input_command_list'

    global sdCoefList ""
    global sdModelIndex 1

    /* loop through commands to get the bootstrapped p-values */
    while "`commands'" != ""{
        
        gettoken cmd commands : commands, parse(" ||") bind match(parenthesis)
                
        qui: store_pvalues `cmd' mat($sdBootPVals)

        /* keep track of the number of models */
        global sdModelIndex = ${sdModelIndex}+1
    }

    matrix rownames $sdBootPVals  = $sdCoefList

    tempname bs_ordered 
    mat_order `bs_ordered' : $sdBootPVals $sdPVal

    /* di "`: rowfullnames `bs_ordered''"   */

    /* stepdown from max to min to enforce monotonicity*/
    forvalues r = `=`n_rows'-1'(-1)1{ // start in the second-to-last row
        /* "cascade" small values down, so that monotonicity is enforced
            - i.e. no entry can be larger than the entries below it */
        mat `bs_ordered'[`r',1] = min(`bs_ordered'[`r'+1,1], `bs_ordered'[`r',1]) 
    }

    mata p_bs_col = st_matrix("`bs_ordered'")
    if `bs_index' == 1 mata p_bs = p_bs_col
    else mata p_bs = p_bs, p_bs_col

 
    /* mat list `bs_ordered' */

    mat drop $sdBootPVals
    restore
}

 

/* mata p_bs */

mata raw_p =  ordered[,1]
mata adj_p =  p_bs :<= raw_p
/* mata adj_p */
mata adj_p = rowsum(adj_p) :/ rownonmissing(adj_p)
/* mata adj_p */

mata st_matrix("$sdPAdjusted", adj_p)

/* step-up */
forvalues r = 2/`n_rows'{ // start in the second row
    /* "cascade" large values up, so that monotonicity is enforced
        - i.e. no entry can be smaller than the entries below it */
    mat $sdPAdjusted[`r',1] = max($sdPAdjusted[`r'-1,1], $sdPAdjusted[`r',1]) 
}

matrix rownames $sdPAdjusted  = $sdCoefList


return matrix p = $sdPVal
return matrix adj_p = $sdPAdjusted

end 

***********************************************************************************************
**************** HELPER FUNCTION TO EXTRACT P-VALUES AFTER RUNNING COMMANDS *******************
***********************************************************************************************

/* internal function */
program store_pvalues 

syntax anything(id="model" name=command everything), ///
    ADJust(string) ///
    mat(name) ///
    [*] // allow for arbitrary command-specific options 
    
    if "`options'" != "" local options , `options'
    `command' `options'

    /* check if the p-value matrix exists */
    cap qui mat list `mat'
    local make_mat = _rc != 0
    
    /* if it doesn't exist, we need to initialize it */
    if `make_mat'{
        local x1 : word 1 of `adjust'

        qui: lincom `x1'
        mat `mat' = (r(p))
    }

     

    foreach var of local adjust{
        global sdCoefList ${sdCoefList} [m${sdModelIndex}:`var']

        if `make_mat' {
            /* skip the first variable if we've just created the matrix using it */
            local make_mat 0
            continue 
        }

        qui lincom `var'
        mat `mat' = (`mat'\r(p))
    }


end 

/* 
FIRST: 
run all the regressions and store the p-values 
+ then sort them so that the matrix is ascending */