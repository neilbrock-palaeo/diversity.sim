



*1. Read in the diversity.sim() fuction and the functions in the legacy\_functions.R file. The code also includes two legacy paleotree functions that have been replaced by more flexible, but much slower versions, so I prefer to use the old ones*





*2. Choose your analysis parameters. Lots to fiddle with here*:

​- local: *number of fossil-bearing localities. Default 200*

\- form: *number of formations (has to be a factor of local, as there will be the same number of localities within each formation. Default 20*

\- form.a: *number of formations which the simulated clade will be allowed to enter. This allows you to simulate for example, that you would not expect to find a marine reptile in a terrestrial formation. Obviously has to be less than form. Default 10*

-ptaph: *the probability that each individual will make it through the taphonomic filter into the fossil record. Can be a number between 0 and 1, or a vector of two numbers between 0 and 1 representing the minimum and maximum probabilities. If the latter,  ptaph will be drawn at random from the range each time bin. Default 0.2*

\- pform: *the probability that each formation makes it through the geological filters. Same options and default as ptaph*

\- ploc: *the probability that each locality will be sampled. Same options  and default as ptaph.*

\- pmist: *the proportion of nodes in the phylogeny estimated incorrectly (phylogeny used in simulating the phylogenetic diversity estimate). Number between 0 and 1. Default 0.1*

\- min.spec: minimum number of species in the simulated clade. Default 1000

\- max.spec: *maximum number of species in the simulated clade. Default 2000*

\- nbins: *number of time bins simulated in the analysis. Default 100*

\- disp: *probability per lineage per time bin of dispersal to a new region. Default 0.25*

\- ext: *probability per lineage per time bin of local extinction from one occupied region. Default 0.25*





*3. For the purposes of this example we'll just use the default settings\*.*



&nbsp;test<-diversity.sim<(local=200,form=20,form.a=10,ptaph=0.2,pform=0.2,ploc=0.2,pmist=0.1,

min.spec=1000,max.spec=2000,nbins=100,disp=0.25,ext=0.25)



test



*\*Note - these parameters are not a recommendation or intended to be realistic; they are chosen purely to provide a reasonably quick illustrative simulation*



*4. The output is a matrix. Each column is a time bin, and each of there are 24 rows:*

\- "True Diversity": *The actual number of species present in each time bin before all sampling filters*

\- "SIB TDE": *Sample in bin taxic diversity estimate. The simplest way to estimate of diversity: count the number of species sampled in each time bin*

\- "RT TDE": *Range through taxic diversity estimate. Again, counting the number of species sampled in each time bin, but also includes Lazarus taxa (unsampled gaps in a species range e.g. a species present in bin 1 and 3 was presumably present in bin 2)*

\- "PDE": *Phylogenetic diversity estimate. Estimating diversity from the phylogeny by counting number of lineages in each time bin (includes ghost lineages; lineages not sampled in the fossil record but inferred from the phylogeny by working on the assumption that two sister taxa must have diverged from their common ancestor at the same time)*

\- "Squares": *The squares equation: S + s1^2 \* sum(N^2) / (sum(N)^2 - s1 \* S), if N is the number of occurrences, S is the number of species, s1 is the number of singletons (species known from one occurrence). Note, some time bins may be Inf; they are the bins where all species were singletons*

\- "RDE Form A": *Residual Diversity Estimate, calculated using number of formations bearing fossils of your clade of interest.*

\- "RDE Form B": *Residual Diversity Estimate, calculated using number of formations bearing fossils of your clade of interest or a larger clade containing your clade of interest (this is an attempt to deal with redundancy: if a clade was in reality more diverse, you'd expect there to be more formations bearing it).*

\- "RDE Form C": *Residual Diversity Estimate, calculated using total number of formations sampled that could have preserved fossils of your clade (i.e. where the clade was permitted to disperse) whether they bear fossils of your clade of interest or not.*

\- "RDE Form D": *Residual Diversity Estimate, calculated using total number of formations sampled whether they bear fossils of your clade of interest or not.*

\- "RDE Loc A": *Residual Diversity Estimate, calculated using number of localities bearing fossils of your clade of interest.*

\- "RDE Loc B": *Residual Diversity Estimate, calculated using number of localities bearing fossils of your clade of interest or a larger clade containing your clade of interest.*

\- "RDE Loc C": *Residual Diversity Estimate, calculated using total number of localities sampled that could have preserved fossils of your clade (i.e. where the clade was permitted to disperse) whether they bear fossils of your clade of interest or not.*

\- "RDE Loc D": *Residual Diversity Estimate, calculated using total number of localities sampled whether they bear fossils of your clade of interest or not.*

\- "Proxy Form A": *Sampling proxy; number of formations bearing fossils of your clade of interest.*

\- "Proxy Form B": *Sampling proxy; number of formations bearing fossils of your clade of interest or a larger clade containing your clade of interest.*

\- "Proxy Form C": *Sampling proxy; total number of formations sampled that could have preserved fossils of your clade (i.e. where the clade was permitted to disperse) whether they bear fossils of your clade of interest or not.*

\- "Proxy Form D": *Sampling proxy;  total number of formations sampled whether they bear fossils of your clade of interest or not.*

\- "Proxy Loc A": *Sampling proxy;  number of localities bearing fossils of your clade of interest.*

\- "Proxy Loc B": *Sampling proxy;  number of localities bearing fossils of your clade of interest or a larger clade containing your clade of interest*

\- "Proxy Loc C": *Sampling proxy;  total number of localities sampled that could have preserved fossils of your clade (i.e. where the clade was permitted to disperse) whether they bear fossils of your clade of interest or not.*

\- "Proxy Loc D": *Sampling proxy;  total number of localities sampled whether they bear fossils of your clade of interest or not.*

\- "pform": *pform for that time bin (will only be variable if the argument was a range)*

\- "ploc": *ploc for that time bin (will only be variable if the argument was a range)*

\- "ptaph": *ptaph for that time bin (will only be variable if the argument was a range)*



