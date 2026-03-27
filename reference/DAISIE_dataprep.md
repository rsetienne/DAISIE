# Prepare colonisation and branching time data to run in DAISIE.

This function produces a data object that can be run in DAISIE
likelihood computation/optimization functions. The function converts a
user-specified table to a DAISIE-compatible format. See
Galapagos_datatable.Rdata for a template of an input table.)

## Usage

``` r
DAISIE_dataprep(
  datatable,
  island_age,
  M,
  number_clade_types = 1,
  list_type2_clades = NA,
  prop_type2_pool = "proportional",
  epss = 1e-05,
  verbose = TRUE
)
```

## Arguments

- datatable:

  Data frame (table) with user-specified data. See file
  Galapagos_datatable.Rdata for a template of an input table. Each row
  on the table represents and independent colonisation event. Table has
  the following four columns.  
    
  `$Clade_name` - name of independent colonization event  
  `$Status` - One of the following categories:  
  \* "Non_endemic": applies to non-endemic species when an approximate
  colonisation time is known  
  \* "Non_endemic_MaxAge": applies to non-endemic species for cases
  where colonisation time is unknown  
  \* "Endemic": applies to endemic species or endemic clades when an
  approximate colonisation time is known  
  \* "Endemic_MaxAge": applies to endemic species or endemic clades for
  cases where the colonisation time is unknown, or when the user wants
  to specify an upper bound for colonisation. This could for example
  apply to endemic species that have recently gone extinct because of
  anthropogenic causes, and which are not included in the phylogeny
  ("NA" should be given in the branching times column). It could also
  apply to insular radiations with long stem branches, for which the
  time of the first cladogenetic event is known, but the precise time of
  colonisation is not.  
  \* "Endemic_MaxAge_MinAge": same as Endemic_MaxAge but also includes a
  minimum age for colonisation.  
  \* "Non_endemic_MaxAge_MinAge": same as Non_endemic_MaxAge but also
  includes a minimum age for colonisation.#'  
  \* "Endemic&Non_Endemic": when endemic clade is present and its
  mainland ancestor has re-colonized  
  `$Missing_species` - Number of island species that were not sampled
  for particular clade (only applicable for "Endemic" clades). If NA is
  given in branching times column, this should be equal to the number of
  species in the clade minus 1  
  `$Branching_times` - Stem age of the population/species in the case of
  "Non_endemic", "Non_endemic_MaxAge" and "Endemic" species with no
  extant close relatives on the island. Set "NA" if colonisation time
  unknown and no upper bound is known. For "Endemic" cladogenetic
  species these should be branching times of the radiation, including
  the stem age of the radiation (colonisation time estimate).  

- island_age:

  Age of island in appropriate units

- M:

  The size of the mainland pool, i.e the number of species that can
  potentially colonize the island

- number_clade_types:

  Number of clade types. Default: number_clade_types = 1 all species are
  considered to belong to same macroevolutionary process. If
  number_clade_types = 2, there are two types of clades with distinct
  macroevolutionary processes.

- list_type2_clades:

  If number_clade_types = 2, list_type2_clades specifies the names of
  the clades that have a distinct macroevolutionary process. The names
  must match those in the \$Clade_name column of the source data table
  (e.g. list_type2_clades = "Finches"). If number_clade_types = 1, then
  list_type2_clades = NA should be specified (default)

- prop_type2_pool:

  Specifies the fraction of potential mainland colonists that have a
  distinct macroevolutionary process. Applies only if number_clade_types
  = 2. Default "proportional" sets the fraction to be proportional to
  the number of clades of distinct macroevolutionary process that have
  colonised the island. Alternatively, the user can specify a value
  between 0 and 1 (e.g. if mainland pool size is 1000 and
  prop_type2_pool = 0.02 then number of type2 species is 20).

- epss:

  Default= 1E-5 should be appropriate in most cases. This value is used
  to set the maximum age of colonisation of "Non_endemic_MaxAge" and
  "Endemic_MaxAge" species to an age that is slightly younger than the
  island for cases when the age provided for that species is older than
  the island. The new maximum age is then used as an upper bound to
  integrate over all possible colonisation times.

- verbose:

  Boolean. States if intermediate results should be printed to console.
  Defaults to `TRUE`.

## Value

- datalist:

  R list object containing data:  
  The first element of the list has two or three components:  
  `$island_age` - the island age  
  Then, depending on whether a distinction between species types is
  made, we have:  
  `$not_present` - the number of mainland lineages that are not present
  on the island  
  or:  
  `$not_present_type1` - the number of mainland lineages of type 1 that
  are not present on the island  
  `$not_present_type2` - the number of mainland lineages of type 2 that
  are not present on the island  
  The following elements of the list each contains information on a
  single colonist lineage on the island and has 5 components:  
  `$colonist_name` - the name of the species or clade that colonized the
  island  
  `$branching_times` - island age and stem age of the population/species
  in the case of "Non-endemic", "Non-endemic_MaxAge" and "Endemic"
  anagenetic species. For "Endemic" cladogenetic species these are
  island age and branching times of the radiation including the stem age
  of the radiation.  
  `$stac` - the status of the colonist  
  \* Non_endemic_MaxAge: 1  
  \* Endemic: 2  
  \* Endemic&Non_Endemic: 3  
  \* Non_endemic: 4  
  \* Endemic_MaxAge: 5 (if only colonisation time was given)  
  \* Endemic_MaxAge: 6 (if colonisation time and cladogenesis times were
  given)  
  `$missing_species` - number of island species that were not sampled
  for particular clade (only applicable for endemic clades)  
  `$type_1or2` - whether the colonist belongs to type 1 or type 2  

## Details

The output is an R list containing the data formatted to be run on other
DAISIE functions.

## References

Valente, L.M., A.B. Phillimore and R.S. Etienne (2015). Equilibrium and
non-equilibrium dynamics simultaneously operate in the Galapagos
islands. Ecology Letters 18: 844-852.

## Author

Luis M Valente

## Examples

``` r



### Create Galapagos data object where all taxa have the same macroevolutionary process

utils::data(Galapagos_datatable, package = "DAISIE")
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000
   )
#> Colonisation time of 7.456 for Coccyzus is older than island age
#> Colonisation time of 10.285 for Pyrocephalus is older than island age
#> [[1]]
#> [[1]]$island_age
#> [1] 4
#> 
#> [[1]]$not_present
#> [1] 992
#> 
#> 
#> [[2]]
#> [[2]]$colonist_name
#> [1] "Coccyzus"
#> 
#> [[2]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[2]]$stac
#> [1] 1
#> 
#> [[2]]$missing_species
#> [1] 0
#> 
#> [[2]]$type1or2
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]]$colonist_name
#> [1] "Dendroica"
#> 
#> [[3]]$branching_times
#> [1] 4.00 0.34
#> 
#> [[3]]$stac
#> [1] 4
#> 
#> [[3]]$missing_species
#> [1] 0
#> 
#> [[3]]$type1or2
#> [1] 1
#> 
#> 
#> [[4]]
#> [[4]]$colonist_name
#> [1] "Finches"
#> 
#> [[4]]$branching_times
#>  [1] 4.0000 3.0282 1.3227 0.8223 0.4286 0.3462 0.2450 0.1180 0.0808 0.0756
#> [11] 0.0527 0.0525 0.0327 0.0322 0.0221 0.0118
#> 
#> [[4]]$stac
#> [1] 2
#> 
#> [[4]]$missing_species
#> [1] 0
#> 
#> [[4]]$type1or2
#> [1] 1
#> 
#> 
#> [[5]]
#> [[5]]$colonist_name
#> [1] "Mimus"
#> 
#> [[5]]$branching_times
#> [1] 4.000 3.958 3.422 2.884 0.459
#> 
#> [[5]]$stac
#> [1] 2
#> 
#> [[5]]$missing_species
#> [1] 0
#> 
#> [[5]]$type1or2
#> [1] 1
#> 
#> 
#> [[6]]
#> [[6]]$colonist_name
#> [1] "Myiarchus"
#> 
#> [[6]]$branching_times
#> [1] 4.000 0.855
#> 
#> [[6]]$stac
#> [1] 2
#> 
#> [[6]]$missing_species
#> [1] 0
#> 
#> [[6]]$type1or2
#> [1] 1
#> 
#> 
#> [[7]]
#> [[7]]$colonist_name
#> [1] "Progne"
#> 
#> [[7]]$branching_times
#> [1] 4.000 0.086
#> 
#> [[7]]$stac
#> [1] 2
#> 
#> [[7]]$missing_species
#> [1] 0
#> 
#> [[7]]$type1or2
#> [1] 1
#> 
#> 
#> [[8]]
#> [[8]]$colonist_name
#> [1] "Pyrocephalus"
#> 
#> [[8]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[8]]$stac
#> [1] 1
#> 
#> [[8]]$missing_species
#> [1] 0
#> 
#> [[8]]$type1or2
#> [1] 1
#> 
#> 
#> [[9]]
#> [[9]]$colonist_name
#> [1] "Zenaida"
#> 
#> [[9]]$branching_times
#> [1] 4.00 3.51
#> 
#> [[9]]$stac
#> [1] 2
#> 
#> [[9]]$missing_species
#> [1] 0
#> 
#> [[9]]$type1or2
#> [1] 1
#> 
#> 

### Create Galapagos data object with a distinct macroevolutionary processes
# for the Darwin's finches. One process applies to type 1 species (all species
# except for Darwin's finches) and the other applies only to type 2 species
# (Darwin's finches). Set fraction of potential colonists of type 2 to be
# proportional to the number of type2 clades present on the island.

utils::data(Galapagos_datatable, package = "DAISIE")
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000,
   number_clade_types = 2,
   list_type2_clades = "Finches"
   )
#> Colonisation time of 7.456 for Coccyzus is older than island age
#> Colonisation time of 10.285 for Pyrocephalus is older than island age
#> [[1]]
#> [[1]]$island_age
#> [1] 4
#> 
#> [[1]]$not_present_type1
#> [1] 868
#> 
#> [[1]]$not_present_type2
#> [1] 124
#> 
#> 
#> [[2]]
#> [[2]]$colonist_name
#> [1] "Coccyzus"
#> 
#> [[2]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[2]]$stac
#> [1] 1
#> 
#> [[2]]$missing_species
#> [1] 0
#> 
#> [[2]]$type1or2
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]]$colonist_name
#> [1] "Dendroica"
#> 
#> [[3]]$branching_times
#> [1] 4.00 0.34
#> 
#> [[3]]$stac
#> [1] 4
#> 
#> [[3]]$missing_species
#> [1] 0
#> 
#> [[3]]$type1or2
#> [1] 1
#> 
#> 
#> [[4]]
#> [[4]]$colonist_name
#> [1] "Finches"
#> 
#> [[4]]$branching_times
#>  [1] 4.0000 3.0282 1.3227 0.8223 0.4286 0.3462 0.2450 0.1180 0.0808 0.0756
#> [11] 0.0527 0.0525 0.0327 0.0322 0.0221 0.0118
#> 
#> [[4]]$stac
#> [1] 2
#> 
#> [[4]]$missing_species
#> [1] 0
#> 
#> [[4]]$type1or2
#> [1] 2
#> 
#> 
#> [[5]]
#> [[5]]$colonist_name
#> [1] "Mimus"
#> 
#> [[5]]$branching_times
#> [1] 4.000 3.958 3.422 2.884 0.459
#> 
#> [[5]]$stac
#> [1] 2
#> 
#> [[5]]$missing_species
#> [1] 0
#> 
#> [[5]]$type1or2
#> [1] 1
#> 
#> 
#> [[6]]
#> [[6]]$colonist_name
#> [1] "Myiarchus"
#> 
#> [[6]]$branching_times
#> [1] 4.000 0.855
#> 
#> [[6]]$stac
#> [1] 2
#> 
#> [[6]]$missing_species
#> [1] 0
#> 
#> [[6]]$type1or2
#> [1] 1
#> 
#> 
#> [[7]]
#> [[7]]$colonist_name
#> [1] "Progne"
#> 
#> [[7]]$branching_times
#> [1] 4.000 0.086
#> 
#> [[7]]$stac
#> [1] 2
#> 
#> [[7]]$missing_species
#> [1] 0
#> 
#> [[7]]$type1or2
#> [1] 1
#> 
#> 
#> [[8]]
#> [[8]]$colonist_name
#> [1] "Pyrocephalus"
#> 
#> [[8]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[8]]$stac
#> [1] 1
#> 
#> [[8]]$missing_species
#> [1] 0
#> 
#> [[8]]$type1or2
#> [1] 1
#> 
#> 
#> [[9]]
#> [[9]]$colonist_name
#> [1] "Zenaida"
#> 
#> [[9]]$branching_times
#> [1] 4.00 3.51
#> 
#> [[9]]$stac
#> [1] 2
#> 
#> [[9]]$missing_species
#> [1] 0
#> 
#> [[9]]$type1or2
#> [1] 1
#> 
#> 

### Create Galapagos data object with a distinct macroevolutionary processes
# for the Darwin's finches. One process applies to type 1 species (all species
# except for Darwin's finches) and the other applies only to type 2 species
# (Darwin's finches). Set fraction of potential colonists of type 2 to be 0.163.

utils::data(Galapagos_datatable, package = "DAISIE")
DAISIE_dataprep(
   datatable = Galapagos_datatable,
   island_age = 4,
   M = 1000,
   number_clade_types = 2,
   list_type2_clades = "Finches",
   prop_type2_pool = 0.163
   )
#> Colonisation time of 7.456 for Coccyzus is older than island age
#> Colonisation time of 10.285 for Pyrocephalus is older than island age
#> [[1]]
#> [[1]]$island_age
#> [1] 4
#> 
#> [[1]]$not_present_type1
#> [1] 830
#> 
#> [[1]]$not_present_type2
#> [1] 162
#> 
#> 
#> [[2]]
#> [[2]]$colonist_name
#> [1] "Coccyzus"
#> 
#> [[2]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[2]]$stac
#> [1] 1
#> 
#> [[2]]$missing_species
#> [1] 0
#> 
#> [[2]]$type1or2
#> [1] 1
#> 
#> 
#> [[3]]
#> [[3]]$colonist_name
#> [1] "Dendroica"
#> 
#> [[3]]$branching_times
#> [1] 4.00 0.34
#> 
#> [[3]]$stac
#> [1] 4
#> 
#> [[3]]$missing_species
#> [1] 0
#> 
#> [[3]]$type1or2
#> [1] 1
#> 
#> 
#> [[4]]
#> [[4]]$colonist_name
#> [1] "Finches"
#> 
#> [[4]]$branching_times
#>  [1] 4.0000 3.0282 1.3227 0.8223 0.4286 0.3462 0.2450 0.1180 0.0808 0.0756
#> [11] 0.0527 0.0525 0.0327 0.0322 0.0221 0.0118
#> 
#> [[4]]$stac
#> [1] 2
#> 
#> [[4]]$missing_species
#> [1] 0
#> 
#> [[4]]$type1or2
#> [1] 2
#> 
#> 
#> [[5]]
#> [[5]]$colonist_name
#> [1] "Mimus"
#> 
#> [[5]]$branching_times
#> [1] 4.000 3.958 3.422 2.884 0.459
#> 
#> [[5]]$stac
#> [1] 2
#> 
#> [[5]]$missing_species
#> [1] 0
#> 
#> [[5]]$type1or2
#> [1] 1
#> 
#> 
#> [[6]]
#> [[6]]$colonist_name
#> [1] "Myiarchus"
#> 
#> [[6]]$branching_times
#> [1] 4.000 0.855
#> 
#> [[6]]$stac
#> [1] 2
#> 
#> [[6]]$missing_species
#> [1] 0
#> 
#> [[6]]$type1or2
#> [1] 1
#> 
#> 
#> [[7]]
#> [[7]]$colonist_name
#> [1] "Progne"
#> 
#> [[7]]$branching_times
#> [1] 4.000 0.086
#> 
#> [[7]]$stac
#> [1] 2
#> 
#> [[7]]$missing_species
#> [1] 0
#> 
#> [[7]]$type1or2
#> [1] 1
#> 
#> 
#> [[8]]
#> [[8]]$colonist_name
#> [1] "Pyrocephalus"
#> 
#> [[8]]$branching_times
#> [1] 4.00000 3.99999
#> 
#> [[8]]$stac
#> [1] 1
#> 
#> [[8]]$missing_species
#> [1] 0
#> 
#> [[8]]$type1or2
#> [1] 1
#> 
#> 
#> [[9]]
#> [[9]]$colonist_name
#> [1] "Zenaida"
#> 
#> [[9]]$branching_times
#> [1] 4.00 3.51
#> 
#> [[9]]$stac
#> [1] 2
#> 
#> [[9]]$missing_species
#> [1] 0
#> 
#> [[9]]$type1or2
#> [1] 1
#> 
#> 
```
