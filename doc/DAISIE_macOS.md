# DAISIE on macOS

This worked on a fresh install of 
* Mojave 10.14.4
* Catalina 10.15.7
* Big Sur 11.2

We install everything via homebrew:

## Install homebrew
```
/bin/bash -c "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/HEAD/install.sh)"
```

## Install gcc and gfortran
```
brew install gcc
```

## Install R
```
brew install r
```

## Add openssl to path
```
echo 'export PATH="/usr/local/opt/openssl@1.1/bin:$PATH"' >> ~/.profile
source .profile
```

## Install libgit2 (devtools needs this)

```
brew install libgit2
```

## Check R version

```
R --version 
R version 4.0.3 (2020-10-10) -- "Bunny-Wunnies Freak Out"
Copyright (C) 2020 The R Foundation for Statistical Computing
Platform: x86_64-apple-darwin18.7.0 (64-bit)
```

## Install RStudio
```
brew install rstudio
```
Optional: install devtools.


## Clean and Rebuild DAISIE.Rproj
```
DONE (DAISIE)
```
