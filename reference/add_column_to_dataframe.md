# Add a column to a data frame

Add a column to a data frame

## Usage

``` r
add_column_to_dataframe(df, position, column_to_insert)
```

## Arguments

- df:

  data frame to add the column to

- position:

  location in data frame where to insert the column. Position can also
  be a name of a column

- column_to_insert:

  the elements of the column to insert. If the column has a name, this
  name will be copied into the data frame. Id is does not have a name,
  it will get the name "nc".

## Value

A data frame with the column inserted
