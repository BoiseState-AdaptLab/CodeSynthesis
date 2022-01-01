# CodeSynthesis
Authors: Tobi Popoola, Flo Ciaglia, Catherine Olschanowsky
## Overview
This is a research piece for generating inspector/executor code for sparse format
conversion. 

## Build Instructions
```shell script
git clone https://github.com/BoiseState-AdaptLab/CodeSynthesis
cd CodeSynthesis
mkdir build/
cd build/
cmake ../
make
```

## Build Scripts
Some preset scripts have been generated by the synthesis algorithm, samples of these conversion 
scripts are available for inspection and check.

``` shell script
   make scripts
```

## Tests
```shell script
cd build/
make test
```

