sudo yum install -y R

sudo yum install git -y

mkdir code
mkdir inputs
mkdir outputs

cd code
git clone https://github.com/jsocolar/colombiaBeta.git     

R

install.packages("remotes")
remotes::install_github("jsocolar/flocker", ref = "more-models")

library(cmdstanr)
install_cmdstan(cores = parallel::detectCores()/2)

cpp_options = list(
  stan_threads=FALSE
  , STAN_NO_RANGE_CHECKS=TRUE
  , CXXFLAGS_OPTIM = "-march=native -mtune=native"
)

cmdstanr::cmdstan_make_local(cpp_options = cpp_options, append = FALSE)
cmdstanr::rebuild_cmdstan(cores = parallel::detectCores()/2)