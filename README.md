# ZILGM
Zero-inflated Local Graphical Models

ZILGM is an R package. ZILGM provides functions to fit the type I and II zero-inflated local negative binomial graphical models (ZILNBGM), and zero-inflated local Poisson graphical model (ZILPGM).
It also provides a function that generates simulation data. 

1. INSTALLATION

The ZILGM package does not list on CRAN repository. Therefore, the ZILGM package can not be installed through install.packages("ZILGM") in R prompt.
The ZILGM package can be installed through our GitHub.
Please install in one of the two installation methods below in R.

(1) From GitHub
```{r}
library(devtools)
install_github("bbeomjin/ZILGM")
```

2. USAGE NOTES

(1) Access to and use of real data in the manuscript

- Cervical data : Cervical data can be available at the "MLSeq" package in R. The data can be loaded with the code below in R.
> library(BioManager)
> BioManager::install("MLSeq")
> library(MLSeq)
> data(cervical)

- Classic3 data : Classic3 data can be available at http://www.dataminingresearch.com/index.php/2010/09/classic3-classic4-datasets. 
	         However, it seems that the URL cannot be accessed. Therefore, we included classic3 data in our ZILGM package.
	         The data can be loaded with the code below in R.
> library(ZILGM)
> data(classic3)


(2) Description of R functions in ZILGM package

- Descriptions of arguments in the functions in ZILGM can be obtained by help() or ? in R prompt, and documentation of ZILGM.   


(3) List of R functions in ZILGM package

- zilgm : "zilgm" function is used to fit the type I and II ZILNBGM model, and ZILPGM model.

- generate_network : "generate_network" function is that generates scale-free, hub and random graph strutures.

- zilgm_sim : "zilgm_sim" is used to generate the zero-inflated count data with overdispersion for simulation given graph structure.

- find_lammax : "find_lammax"  function finds the maximum value of the regularization parameter.


3. EXAMPLE

# Generation of simulated data under the random graph structure
> require(ZILGM)
> set.seed(1)
> n = 100; p = 10; prob = 2 / p;
> A = generate_network(p, prob, type = "random")
> simul_dat = zilgm_sim(A = A, n = n, p = p, zlvs = 0.1, family = "negbin", signal = 1.5, theta = 0.25, noise = 0.0)    

# Compute a sequence of regularization parameter
> lambda_max = find_lammax(simul_dat$X)
> lambda_min = 1e-4 * lambda_max
> lambs = exp(seq(log(lambda_max), log(lambda_min), length.out = 50))
> nb2_fit = zilgm(X = simul_dat$X, lambda = lambs, family = "NBII", update_type = "IRLS", do_boot = TRUE,
                      boot_num = 30, sym = "OR")

# Get estimated graph
> est_graph = nb2_fit$network[[nb2_fit$opt_index]]
