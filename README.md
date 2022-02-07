# NIPET/CuVec speed tests

Testing https://github.com/NiftyPET/NIPET OSEM reconstruction before and after https://github.com/AMYPAD/CuVec integration.

1. Create two different Python3 environments
    - e.g. using `conda create -n old python=3.9` and `conda create -n new python=3.9`
2. Configure paths in the `Makefile`
    - `PATHTOOLS`: Any directory in which to store extra binaries
    - `HMUDIR`: Directory containing proprietary Siemens hardware attenuation maps
    - `DATA`: Directory containing PET scan raw (listmode) data (e.g. from <https://doi.org/10.5281/zenodo.3877529>)
    - `old.%: PY`: Path to Python3 executable (for testing old version)
    - `new.%: PY`: Path to Python3 executable (for testing new version)
3. Setup
    - `make old.setup`
    - `make new.setup`
4. Run
    - `make old.cprof`
    - `make new.cprof`
