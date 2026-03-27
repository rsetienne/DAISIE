# Calculates information criterion from DAISIE ML estimates?

Calculates information criterion from DAISIE ML estimates?

## Usage

``` r
DAISIE_IC(
  datalist,
  initparsopt,
  idparsopt,
  parsfix,
  idparsfix,
  endmc = 1000,
  res = 100,
  cond = 0,
  ddmodel = 0
)
```

## Arguments

- datalist:

  Data object containing information on colonisation and branching
  times. This object can be generated using the DAISIE_dataprep
  function, which converts a user-specified data table into a data
  object, but the object can of course also be entered directly. It is
  an R list object with the following elements.  
  The first element of the list has two or three components:  
    
  `$island_age` - the island age  
  Then, depending on whether a distinction between types is made, we
  have:  
  `$not_present` - the number of mainland lineages that are not present
  on the island  
  or:  
  `$not_present_type1` - the number of mainland lineages of type 1 that
  are not present on the island  
  `$not_present_type2` - the number of mainland lineages of type 2 that
  are not present on the island  
    
  The remaining elements of the list each contains information on a
  single colonist lineage on the island and has 5 components:  
    
  `$colonist_name` - the name of the species or clade that colonized the
  island  
  `$branching_times` - island age followed by stem age of the
  population/species in the case of Non-endemic, Non-endemic_MaxAge
  species and Endemic species with no close relatives on the island. For
  endemic clades with more than one species on the island (cladogenetic
  clades/ radiations) these should be island age followed by the
  branching times of the island clade including the stem age of the
  clade  
  `$stac` - the status of the colonist  
    
  - Non_endemic_MaxAge: 1  
  - Endemic: 2  
  - Endemic&Non_Endemic: 3  
  - Non_Endemic: 4  
  - Endemic_Singleton_MaxAge: 5  
  - Endemic_Clade_MaxAge: 6  
  - Endemic&Non_Endemic_Clade_MaxAge: 7  
  - Non_endemic_MaxAge_MinAge: 8  
  - Endemic_Singleton_MaxAge_MinAge: 9  
    
  `$missing_species` - number of island species that were not sampled
  for particular clade (only applicable for endemic clades)  
  `$type1or2` - whether the colonist belongs to type 1 or type 2  

- initparsopt:

  The initial values of the parameters that must be optimized, they are
  all positive.

- idparsopt:

  The ids of the parameters that must be optimized. The ids are defined
  as follows:  
    
  id = 1 corresponds to lambda^c (cladogenesis rate)  
  id = 2 corresponds to mu (extinction rate)  
  id = 3 corresponds to K (clade-level carrying capacity)  
  id = 4 corresponds to gamma (immigration rate)  
  id = 5 corresponds to lambda^a (anagenesis rate)  
  id = 6 corresponds to lambda^c (cladogenesis rate) for an optional
  subset of the species  
  id = 7 corresponds to mu (extinction rate) for an optional subset of
  the species  
  id = 8 corresponds to K (clade-level carrying capacity) for an
  optional subset of the species  
  id = 9 corresponds to gamma (immigration rate) for an optional subset
  of the species  
  id = 10 corresponds to lambda^a (anagenesis rate) for an optional
  subset of the species  
  id = 11 corresponds to p_f (fraction of mainland species that belongs
  to the second subset of species.

- parsfix:

  The values of the parameters that should not be optimized.

- idparsfix:

  The ids of the parameters that should not be optimized, e.g. c(1,3) if
  lambda^c and K should not be optimized.

- endmc:

  Numeric for how many simulations should run.

- res:

  Sets the maximum number of species for which a probability must be
  computed, must be larger than the size of the largest clade.

- cond:

  cond = 0 : conditioning on island age  
  cond = 1 : conditioning on island age and non-extinction of the island
  biota  
  . cond \> 1 : conditioning on island age and having at least cond
  colonizations on the island. This last option is not yet available for
  the IW model  

- ddmodel:

  Sets the model of diversity-dependence:  
    

  - ddmodel = 0 : no diversity dependence

  - ddmodel = 1 : linear dependence in speciation rate

  - ddmodel = 11: linear dependence in speciation rate and in
    immigration rate

  - ddmodel = 2 : exponential dependence in speciation rate

  - ddmodel = 21: exponential dependence in speciation rate and in
    immigration rate

## Value

List of two numerics with WIC and AICb
