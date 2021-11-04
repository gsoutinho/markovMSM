# markovMSM: An R package for checking the Markov condition in multi-state survival data

  ```markovMSM``` is an R package which considers tests of the Markov assumption 
that are applicable to general multi-state models. Three approaches using existing
methodology are considered: a simple method based on including
covariates depending on the history in Cox models for the transition intensities; methods
based on measuring the discrepancy of the non-Markov estimators of the transition probabilities
to the Markovian Aalen-Johansen estimators; and, finally, methods that were
developed by considering summaries from families of log-rank statistics where patients
are grouped by the state occupied of the process at a particular time point.

```Installation```If you want to use the release version of the markovMSM package, you can install the package from CRAN as follows:
install.packages(pkgs="markovMSM");

```Authors``` Gustavo Soutinho and Luís Meira-Machado lmachado@math.uminho.pt 
Maintainer: Gustavo Soutinho gustavosoutinho@sapo.pt

```Funding``` This research was financed by Portuguese Funds through FCT - “Fundação para a Ciência e a
Tecnologia", within Projects projects UIDB/00013/2020, UIDP/00013/2020 and the research
grant PD/BD/142887/2018.

```References``` 
Aalen O, Johansen S (1978). “An Empirical transition matrix for non homogeneous Markov
and chains based on censored observations.” Scandinavian Journal of Statistics, 5, 141–150.

Andersen P, Esbjerg S, Sorensen T (2000). “Multistate models for bleeding episodes and
mortality in liver cirrhosis.” Statistics in Medicine, (19), 587–599.

Andersen P, Keiding N (2002). “Multi-state models for event history analysis.” Statistical
Methods in Medical Research, (11), 91–115.

Andersen PK, Borgan Ø, Gill RD, Keiding N (1993). Statistical Models Based on Counting
Processes. Springer-Verlag, New York.

Borgan O (2005). Encyclopedia of biostatistics: Aalen-Johansen estimator. John Wiley &
Sons.

Chiou S, Qian J, Mormino E, Betensky R (2018). “Permutation tests for general dependent
truncation.” Computational Statistics & Data Analysis, 318, 308–324. doi:10.1016/j.
csda.2018.07.012.

Datta S, Satten G (2001). “Validity of the Aalen-Johansen estimators of stage occupation
probabilities and Nelson Aalen integrated transition hazards for non-Markov models.”
Statistics & Probability Letters, 55, 403–411.

de Uña-Álvarez J, Meira-Machado L (2015). “Nonparametric estimation of transition probabilities
in the non-Markov illness-death model: A comparative study.” Biometrics, 71(2),
364–375. ISSN 0006-341X.

Hougaard P (2000). Analysis of Multivariate Survival Data. Statistics for Biology and Health.
Springer-Verlag, New York.

Kay R (1986). “A Markov model for analyzing cancer markers and disease states in survival
studies.” Biometrics, (42), 457–481.
Meira-Machado L, de Uña-Álvarez J, Cadarso-Suárez C (2006). “Nonparametric Estimation
of Transition Probabilities in a Non-Markov Illness-Death Model.” Lifetime Data Analysis,
12, 325–344.