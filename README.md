# swnsGreySeal - Bayesian Population Model for Southwest Nova Scotia Grey Seals

This repository contains a Shiny application for fitting and tuning a Bayesian population model for Southwest Nova Scotia (SWNS) Grey Seals in Template Model Builder (TMB). This work is based on a previous analysis by Steven Rossi, Yanjun Wang, Nell den Heyer, and Hugues Benôit.

## Usage

First, clone the repository to your local machine:

```R
git clone https://github.com/stevenrossi/swnsGreySeal
```

Three packages are needed to run the Shiny app (`shiny`,`TMB`,`tmbstan`), so make sure those are installed:

```R
install.packages( c("shiny","TMB","tmbstan") )
```

Then navigate to the `swnsGreySeal` directory in R and source `server.R` to load the necessary packages and compile and load the TMB model. Model compilation may take a minute or two. Once the model is loaded, we can launch the Shiny app:

```R
shiny::runApp()
```

Upon launching, the model will be fit under default settings, which should take less than 20 seconds. Any of the parameters, priors, MCMC controls, or projection settings on the left-hand side can then be changed, which will immediately trigger a model re-fit. MCMC progress can be tracked in the R console.

