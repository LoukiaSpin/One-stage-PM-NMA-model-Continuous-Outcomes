---
output:
  github_document
---



```{r setup, include = FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

***

#  Modelling observed and missing continuous outcome data in one-stage pattern-mixture network meta-analysis

## The pattern-mixture model

Mavridis et al. [1] first proposed a two-stage approach that, first, models the observed and missing continuous outcomes simultaneously in each arm of every trial via a pattern-mixture model to obtain adjusted within-trial results (treatment effect and standard error), and then, it pools the trials using a pairwise or network meta-analysis model. The pattern-mixture model distinguishes the participants to those who remained to the trial and those who left the trial early for any reason. Then, in each subgroup, the mean and variance of the outcome is considered. However, this information is not known in the subgroup with the missing participants (except if these participants have been followed-up after leaving the trial).

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


## Installation

We have developed the R function <span style="color: blue;">`NMA.IMDOM.IMROM.PM.model()`</span> to handle MCOD in network meta-analysis efficiently using the informative missingness parameters and effect measures described above. The IMDoM parameter is intuitively related to the mean difference and the standardised mean difference, whereas the IMRoM parameter can be used in conjunction with the ratio of means in the logarithmic scale. Before using the R function, it is necessary to install the libraries `dplyr` and `R2jags` to allow the function to perform a required data management and to implement Bayesian analysis in [JAGS](https://sourceforge.net/projects/mcmc-jags/) (in case JAGS is not downloaded yet). 

You can download the __package__ with the proposed function directly from [GitHub](https://github.com/LoukiaSpin/One-stage-PM-NMA-model-Continuous-Outcomes/) or use the R package [devtools](https://cran.r-project.org/web/packages/devtools/index.html):

```r
install.packages("devtools")
devtools::install_github("LoukiaSpin/Quantifying-Robustness-in-Meta-analysis", build_vignettes = T)
```
To obtain the necessary results for all assumptions and structures about the missingness parameter, we have intergrated the function `collect jags results Full NMA PM Continuous_function()` in our main function `NMA.IMDOM.IMROM.PM.model()`.

### Required format of the dataset 

To use this R function, the dataset must have a wide-format structure, where every trial occupies one row and every intervention-arm occupies one column. An additional column, *NMA*, is needed to indicate the number of networks that comprise the dataset - even if there is only one network in the dataset.


Then, the output is a series of text-files and it is saved __automatically__ in the folder <span style="color: blue;">__Output__</span>:

* Pattern-mixture IMDOM & IMROM model
  + Model scripts
    + Full NMA model
      + Fixed-effect model
        + Mean Difference
        + Ratio of Means
        + Standardised Mean Difference
      + Random-effects model
        + Mean Difference
        + Ratio of Means
        + Standardised Mean Difference
  + <span style="color: blue;">__Output__</span>
    + <span style="color: blue;">Data</span>
      + <span style="color: blue;">Full NMA models</span>
        + <span style="color: blue;">Fixed-effect model</span>
          + <span style="color: blue;">Mean Difference</span>
          + <span style="color: blue;">Ratio of Means</span>
          + <span style="color: blue;">Standardised Mean Difference</span>
        + <span style="color: blue;">Random-effects model</span>
          + <span style="color: blue;">Mean Difference</span>
          + <span style="color: blue;">Ratio of Means</span>
          + <span style="color: blue;">Standardised Mean Difference</span>
  + R scripts
    + Functions
      + Full NMA model


The function creates automatically a sub-folder for each assumption, where it 'drops' the corrresponding text files that refer to the results on the pooled mean difference for all possible comparisons, the within-trial mean difference, the common between-trial variance, the SUCRA values of the interventions, the order of the interventions, and the estimated missingness parameter according to the assumption. For the aforementioned parameters, we obtain the posterior distribution as provided by the `jags()` function alongside the Rubin and Gelman Rhat statistics. 

Currently, the R function `NMA.IMDOM.IMROM.PM.model()` displays a list of results on the aforementioned model parameters for each assumption about the missingness parameter. 
We plan to replace this output with proper illustration, such as forestplots for the pooled effect sizes and barplots for the SUCRA values to visualise the results under all assumptions simultaneously.


## References
1. Mavridis D, White IR, Higgins JP, Cipriani A, Salanti G. Allowing for uncertainty due to missing continuous outcome data in pairwise and network meta-analysis. Stat Med. 2015;34(5):721-41. doi: 10.1002/sim.6365.
2. Spineli LM. An empirical comparison of Bayesian modelling strategies for missing binary outcome data in network meta-analysis. BMC Med Res Methodol. 2019;19(1):86. doi: 10.1186/s12874-019-0731-y.
3. Turner NL, Dias S, Ades AE, Welton NJ. A Bayesian framework to account for uncertainty due to missing binary outcome data in pairwise meta-analysis. Stat Med. 2015;34(12):2062-80.
4. Rhodes KM, Turner RM, Higgins JP. Predictive distributions were developed for the extent of heterogeneity in meta-analyses of continuous outcome data. J Clin Epidemiol. 2015;68(1):52-60. doi: 10.1016/j.jclinepi.2014.08.012.
5. Nikolakopoulou A, Chaimani A, Veroniki AA, Vasiliadis HS, Schmid CH, Salanti G. Characteristics of networks of interventions: a description of a database of 186 published networks. PLoS One. 2014;9(1):e86754. doi: 10.1371/journal.pone.0086754.
6. Friedrich JO, Adhikari NK, Beyene J. Ratio of geometric means to analyze continuous outcomes in meta-analysis: comparison to mean differences and ratio of arithmetic means using empiric data and simulation. Stat Med. 2012;31(17):1857-86. doi: 10.1002/sim.4501.
