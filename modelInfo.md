# Southwest Nova Scotia Grey Seal population model

## Background

Grey seals [(*Halichoerus grypus*)](https://en.wikipedia.org/wiki/Grey_seal) in the northwest Atlantic Ocean comprise a single population with major breeding colonies at Sable Island and in the Gulf of St. Lawrence, and smaller colonies occurring along the coasts of Nova Scotia and the northeastern United States. Inferring abundance at the smaller grey seal colonies is difficult as observations at these colonies have greater uncertainty, comprise relatively short time series and dynamics are presumed to be driven by emigration from larger adjacent colonies. Recent assessment models for grey seals in Canadian waters have aggregated Sable Island and coastal Nova Scotia seals into a single management unit, thereby avoiding direct inferences about dynamics at the smaller colonies.

Atlantic cod [(*Gadus morhua*)](https://en.wikipedia.org/wiki/Atlantic_cod) stocks in the Northwest Atlantic declined to low abundance in the late 1980s / early 1990s due to overfishing and have largely failed to recover despite large scale reductions in fishing effort. Grey seal predation on adult cod appears to be the primary impedement to cod recovery in the southern Gulf of St. Lawrence. Links between grey seal predation and elevated cod mortality have been hypothesized in other ecosystems, such as the Scotian Shelf and Georges Bank. Evaluating these hypotheses requires finer-scale grey seal abundance estimates than are currently available from existing assessment models. To this end, we developed a Bayesian population model to infer grey seal abundance at colonies along southwest Nova Scotia (SWNS).

## Data

Pop counts (*n*) at SWNS colonies (Mud Is., Round Is., Noddy Is., and Flat Is.,) were available for 2007 (*n = 204*), 2010 (*n=417*), 2016 (*n=1849*), and 2021 (*n=2246*). The first two observations were visual counts in late January of mostly weaned pups, with many pups thought to have dispersed given the few mothers left at the colonies.

Raw pup counts (*n*) were converted to total pup production (*I*) by accounting for uncounted pups.
For our analysis, we inflated the 2007 and 2010 counts by 50% to account for pups that were missed in the visual counts or had dispersed before the survey. The latter two observations were based on aerial photographic surveys in early January, which were corrected for the small proportion of seals born after the survey date.

Coefficients of variation (CVs) from statistical analyses were only available for the 2017 and 2021 observations. We assigned large CVs (0.50) to the 2007 and 2010 observations as these counts are highly uncertain. 

| Year | Raw count (*n*) &emsp;	| Adjusted count (*I*) &emsp; | CV  &emsp; 	|
|------	|------	|----------	|------	|
| 2007  &emsp;&emsp;&emsp;  | 204  &emsp;&emsp;&emsp;	| 306   &emsp;&emsp;&emsp;   	| 0.50 	|
| 2010  &emsp;&emsp;&emsp;  | 417  &emsp;&emsp;&emsp;	| 626   &emsp;&emsp;&emsp;   	| 0.50 	|
| 2016  &emsp;&emsp;&emsp;  | 1849 &emsp;&emsp;&emsp;	| 2105  &emsp;&emsp;&emsp;   	| 0.07 	|
| 2021  &emsp;&emsp;&emsp;  | 2246 &emsp;&emsp;&emsp;	| 2420  &emsp;&emsp;&emsp;   	| 0.08 	|

## Population model

Annual pup production at SWNS colonies is represented by
a logistic population model, i.e.,

`I[t] = K*I[0]/( I[0] + (K-I[0])*exp(-r*t) )`

with the following notation:

- `I` - Pup production
- `K` - Carrying capacity (estimated)
- `r` - Population growth rate (estimated)
- `t` - Time step (years)

Given the sparsity of the data, we were unable to estimate all three model parameters (`r`, `K`, and `I[0]`), so we fixed `I[0]` at 10, which was chosen to represent a small initial value

We note that this model assumes that the population is closed to immigration or emigration and that `r` represents an intrinsic population growth rate. These assumptions are presumably violated, given that the growth of SWNS is thought to be driven by emigration from Sable Island. However, estimating a emigration rate from Sable Island to SWNS in not feasible due to the lack of data, so the logistic model at least offers a tractable first step for estimating abundance.

## Objective function

### Likelihood

Observed pup counts in year `t` (`Iobs[t]`) were assumed to arise from log-normal distributions:

`log(Iobs[t]) ~ N( log(I[t]), Isd[t]^2 )`

where variance `Isd^2` was set to `(CV[t]^2 + 1)`.

### Priors

Log-normal priors are specified for `r` and `K`:

`log(r) ~ N( log(rMu), rSD )`

`log(K) ~ N( log(KMu), KSD )`

`rMu`, `rSD`, `KMu`, and `KSD` can all be specified by the user. By default, we set `rMu` to 0.2 and `rSD` to 0.5, which produces a prior with a mode of `exp(log(0.2)-0.5^2) = 0.156` and a 95th percentile of 0.456. We note that a sustained growth rate of 12% was estimated for grey seals on Sable Island, which could be considered an upper bound for the species given the extremely favourable conditions experienced by grey seals during this period. However, given that the growth rate in our model is driven by immigration, it is likely that the rate of growth at SWNS colonies has exceeded 12%.

Priors for `K` are more difficult to determine. Ideally we want to use a wide prior that has minimal influence. We chose default values of `kMu` = 5 and `kSD` = 10.

## Implementation

The statistical model described above was implemented using the Template Model Builder (`tmb`) package within R. Bayes posterior distributions for parameters and predictive pup count distributions were generated using [Hamiltonian Monte Carlo (HMC)](https://arxiv.org/abs/1701.02434), which approximates the posterior density by simulating the evolution of a Hamiltonian system. A key feature of HMC is its ability to move between distant regions of the target density in a single transition using the gradient of the density, thus avoiding the inefficient random walk mechanics of more traditional MCMC algorithms such as Metropolis-Hastings. The algorithm can be [visualized](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12681) as a frictionless disk gliding over a surface, where the position of the disk represents a set of parameter values and the height or potential energy of the disk is analogous to the negative-log-posterior density. The sum of the kinetic and potential energy of the disk, known as the Hamiltonian (H), should be constant as the surface is frictionless, however, paths must be approximated via numerical techniques, causing H to vary. Excessive posterior curvature or transition size can cause the simulated H to diverge from the true H. The accuracy of the HMC algorithm is not guaranteed in the presence of these so-called “divergent transitions”, so ensuring proper posterior geometry and tuning the transition size are important steps in generating reliable inference.


We used an HMC algorithm known as the [no-U-turn sampler (NUTS)](https://jmlr.org/papers/volume15/hoffman14a/hoffman14a.pdf), which eliminates a manual tuning stage in the original HMC algorithm while retaining or improving efficiency. NUTS is accessed via the `tmbstan` R package.

Convergence is monitored using (i) the potential scale reduction factor on rank-normalized split chains (Rhat), where Rhat > 1.01 indicates problems with convergence, (ii) the effective sample size (ESS) of the rank-normalized draws, where an ESS of at least 300 (100 per chain) was considered acceptable, and (iii) the minimum energy Bayesian fraction of missing information (BFMI), which should be at least 0.3 for each chain.
 
[HMC-specific and NUTS-specific diagnostics]((https://mc-stan.org/misc/warnings.html)) (number of divergent transitions of number of transitions that saturated the maximum tree depth, respectively) are also monitored. 

The user specifies the number of chains (`nChain`), number of iterations (`nIter`), and two NUTS settings: `adapt_delta` and `max_treedepth`. The first half of each chain is discarded as a warm-up, so the final number of posterior samples is `nIter*nChain/2`.


## Scaling to total abundance

Total SWNS grey abundance is calculated by scaling estimated pup production according to the user-specified pup-to-adult ratio, which was calculated as 0.26 for the Sable Island population.


## Projection

While the model is only fitted to data up to 2021, the model can be projected for an arbitrary number of years by adjusting the `Final year` control. The vertical line on the time-series plots indicates the beginning of the projection period.



## References

den Heyer, C.E., Lang, S.L.C., Bowen, W.D., and Hammill, M.O. 2017. Pup Production at Scotian Shelf Grey Seal (Halichoerus grypus) Colonies in 2016. DFO Can. Sci. Advis. Sec. Res. Doc. 2017/056. v + 34 p.

den Heyer, C, A. Mosnier, GB Stenson, D. Lidgard, W.D. Bowen and M.O. Hammill . 2021. Pup production of Northwest Atlantic grey seals in Canadian waters. DFO Can. Sci. Advis. Sec. Res Doc. 2021/

Duane, S., Kennedy, A.D., Pendleton, B.J., and Roweth, D. 1987. Hybrid monte carlo. Phys. Lett. 195(2): 216–222. doi:10.1016/0370-2693(87)91197-X.

Hammill, M.O., Stenson, G.B., Proust, F., Carter, P., McKinnon, D., 2007. Feeding by grey seals in the Gulf of St. Lawrence and around Newfoundland. NAMMCO Scient. Publ. 6, 135–152.

Hammill, M.O., Stenson, G.B., Swain, D.P., Benoît, H.P., 2014. Feeding by grey seals on endangered stocks of Atlantic cod and white hake. ICES J. Mar. Sci. 71, 1332–1341. doi:10.1093/icesjms/fsu123.

Hammill, M.O., Heyer, den, C.E., Bowen, W.D., Lang, S.L.C. 2017. Grey Seal Population Trends in Canadian Waters, 1960-2016 and Harvest Advice. DFO Canadian Science Advisory Secretariat Research Document 2017/052. v + 30 p. 

Hoffman, M.D., and Gelman, A. 2014. The No-U-Turn Sampler: Adaptively Setting Path Lengths in Hamiltonian Monte Carlo. J Mach Learn Res (15): 1593–1623.

Kristensen, K. 2018. tmbstan: MCMC Sampling from 'TMB' Model Object using ‘Stan’. R package version 1.0.1. https://CRAN.R-project.org/package=tmbstan.

Kristensen, K., Nielsen, A., Berg, C.W., Skaug, H., and Bell, B.M. 2016. TMB: Automatic Differentiation and Laplace Approximation. J. Stat. Soft. 70(5): 1–21. doi:10.18637/jss.v070.i05.

Monnahan, C.C., Thorson, J.T., Branch, T.A. 2017. Faster estimation of Bayesian models in ecology using Hamiltonian Monte Carlo. Methods in Ecology and Evolution. 8(3):339-348.

Neal, R.M. 2011. MCMC Using Hamiltonian Dynamics. In Handbook of Markov Chain Monte Carlo. pp. 113–162. Ed. by Brooks, S., Gelman, A., Jones G.L. & Meng, X.-L.)

R Core Team. 2022. R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. https://www.R-project.org/.

Rossi, S.P., Cox, S.P., Hammill, M.O., den Heyer, C.E. , Swain, D.P., Mosnier, A., Benoît, H.P. 2021. Forecasting the response of a recovered pinniped population to sustainable harvest strategies that reduce their impact as predators. ICES Journal of Marine Science 78:1804–1814. doi:10.1093/icesjms/fsab088

Vehtari, A., Gelman, A., Simpson, D., Carpenter, B., & Bürkner, P. C. 2021. Rank-normalization, folding, and localization: An improved R ̂ for assessing convergence of MCMC (with discussion). Bayesian analysis, 16(2), 667-718.










