***

#  Modelling observed and missing continuous outcome data in one-stage pattern-mixture network meta-analysis

## The pattern-mixture model

Mavridis et al. [1] have proposed a two-stage approach that, first, models the observed and missing continuous outcomes simultaneously in each arm of every trial via a pattern-mixture model to obtain adjusted within-trial results (treatment effect and standard error), and then, it pools the trials using a pairwise or network meta-analysis model. The pattern-mixture model distinguishes the participants to those who remained to the trial and those who left the trial early for any reason. Then, in each subgroup, the mean and variance of the outcome is considered. However, this information is not known in the subgroup with the missing participants (except if these participants have been followed-up after leaving the trial).

### Two missingness parameters

Mavridis et al. [1] proposed two missingness parameters for MCOD that faciliate the determination of an informative prior distribution that reflects our prior beliefs about the missingness mechanism. 

#### <ins>Informative missingness difference of means (IMDoM)</ins>

The IMDoM parameter is defined as the difference between the mean outcome among missing participants and the mean outcome among those completing the trial

#### <ins>Informative missingness ratio of means (IMRoM)</ins>

The IMRoM parameter is defined as the ratio of the mean outcome among the missing participants to the mean outcome among those completing the trial.

Then, a normal prior distribution is a natural choice for both missingness parameters, where the mean implies our belief about the missingness process on average, and the variance indicates our uncertainty about this prior belief. Specifically:

+ a positive mean indicates that a larger outcome on average is more likely to occur among missing participants rather than completers;
+ a negative mean indicates that a larger outcome on average is less likely to occur among missing participants rather than completers; and
+ a zero mean indicates the missing at random (MAR) assumption on average.


### Structural assumptions for the missingness parameters

The following three assumptions can be made for the missingness parameters [2,3]:

+ *common-within-network* assumption implies the missingness mechanism to be the same in the whole network, and only one parameter is estimated per network;
+ *trial-specific* assumption implies that the missingness mechanisms are different across trials but assumed to be the same in the compared arms resulting in down-weighting trials with unbalanced MCOD. This assumption estimates as many missingness parameters as the number of included trials;
+ *intervention-specific* assumption allows the missingness parameter to be different across interventions but shared across trials, thus resulting in down-weighting trials with higher total MCOD. This assumption estimates as many missingness parameters as the number of interventions in a network.

Then, for each assumption, the prior distribution of the missingness parameter can be structured to be either identical (i.e. constant) or hierarchical (i.e. exhangeable). An additional, more flexible structure is the *independent structure* where the missingness mechanism is considered to differ across arms and trials. This corresponds to estimating one missingness parameter for each arm of every trial. Then, under the independent structure, the missingness parameters may be further assumed to be either *uncorrelated* or *correlated* across the arms of every trial with correlation parameter equal to 0.5 in each trial-arm (following Mavridis et al. [1]).   

## Effect measures for continuous outcome 

In systematic reviews, mean difference and standardised mean difference are the most popular effect measures for continuous outcomes [4, 5]. A less popular effect measure for continuous outcomes is the ratio of (arithmetic) means which has been advocated for its ease of interpretation (as compared to  standardised mean difference) and its ability to pool trials with outcomes in different units (contrary to mean difference) [6]. 


## Description of the repository

The repository includes separate folders for data, models and R (code/scripts), respectively.
* The _data_ folder includes two input text files: 21370258_Stowe(2011).txt and 24996616_Schwingshackl(2014).txt. These are the examples that we considered in our article.
* The _model_ folder includes text files on the proposed Bayesian random-effects network meta-analysis models for aggregate outcomes with missing outcome data;
* The _R_ folder includes the function <span style="color: blue;">`run.model()`</span> to handle MCOD in network meta-analysis efficiently using the informative missingness parameters and effect measures described above. It also includes the script <span style="color: blue;">`Function application`</span> to run both examples using the function<span style="color: blue;">`run.model()`</span> straightfowardly. The IMDoM parameter is intuitively related to the mean difference and the standardised mean difference, whereas the IMRoM parameter can be used in conjunction with the ratio of means in the logarithmic scale. 

Before using the R function, it is necessary to install the libraries `dplyr` and `R2jags` to allow the function to perform a required data management and to implement Bayesian analysis in [JAGS](https://sourceforge.net/projects/mcmc-jags/) (in case JAGS is not downloaded yet). 


### Perform Bayesian random-effects network meta-analysis to handle aggregate MCOD

```r
run.model(data, measure, assumption, mean.misspar, var.misspar, D, n.chains, n.iter, n.burnin, n.thin)
```

#### Explaining the arguments

* data: A data-frame of a one-trial-per-row format containing arm-level data of each trial. This format is widely used for BUGS models. 
* measure: Character string indicating the effect measure with values "MD", "SMD", or "ROM".
* assumption: Character string indicating the structure of the informative missingness parameter. Set assumption equal to one of the following: "HIE-COMMON", "HIE-TRIAL", "HIE-ARM", "IDE-COMMON", "IDE-TRIAL", "IDE-ARM", "IND-CORR", or "IND-UNCORR".
* mean.misspar: A positive non-zero number for the mean of the normal distribution of the informative missingness parameter.
* var.misspar: A positive non-zero number for the variance of the normal distribution of the informative missingness parameter.
* D: A binary number for the direction of the outcome. Set D = 1 for a positive outcome and D = 0 for a negative outcome.
* n.chains: Integer specifying the number of chains for the MCMC sampling; an argument of the jags function in R2jags.
* n.iter: Integer specifying the number of Markov chains for the MCMC sampling; an argument of the jags function in R2jags.
* n.burnin: Integer specifying the number of iterations to discard at the beginning of the MCMC sampling; an argument of the jags function in R2jags.
* n.thin: Integer specifying the thinning rate for the MCMC sampling; an argument of the jags function in R2jags.


#### Output of the function

The function returns results on the pooled effect estimate (mean difference, standardised mean difference, or ratio of means in the log scale) for all possible comparisons, the within-trial effect estimate, the common between-trial standard deviation, the SUCRA values of each intervention, the ranking probability of each intervention for every ranking, and the estimated missingness parameter according to the assumption. Furthermore, it provides results for the predictions of all possible comparisons. For the aforementioned parameters, we obtain the posterior distribution as provided by the `jags()` function alongside the Rubin and Gelman Rhat statistics. 

Currently, the R function `run.model()` displays a list of results on the aforementioned model parameters for each assumption about the missingness parameter. 
We plan to replace this output with proper illustration, such as forestplots for the pooled effect sizes and barplots for the SUCRA values.


## References
1. Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for uncertainty due to missing continuous outcome data in pairwise and network meta-analysis. Stat Med. 2015;34(5):721-41. doi: 10.1002/sim.6365.
2. Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. BMC Med Res Methodol. 2019;19(1):86. doi: 10.1186/s12874-019-0731-y.
3. Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for uncertainty due to missing binary outcome data in pairwise meta-analysis. Stat Med. 2015;34(12):2062-80.
4. Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed for the extent of heterogeneity in meta-analyses of continuous outcome data. J Clin Epidemiol. 2015;68(1):52-60. doi: 10.1016/j.jclinepi.2014.08.012.
5. Nikolakopoulou A, Chaimani A, Veroniki AA, Vasiliadis HS, Schmid CH, Salanti G. Characteristics of networks of interventions: a description of a database of 186 published networks. PLoS One. 2014;9(1):e86754. doi: 10.1371/journal.pone.0086754.
6. Friedrich JO, Adhikari NK, Beyene J. Ratio of geometric means to analyze continuous outcomes in meta-analysis: comparison to mean differences and ratio of arithmetic means using empiric data and simulation. Stat Med. 2012;31(17):1857-86. doi: 10.1002/sim.4501.
