# Package Repository | _SAEN-LARS_    [![DOI](https://zenodo.org/badge/125818623.svg)](https://zenodo.org/badge/latestdoi/125818623)

* __Version__: 1.0.0 
* __Title__: Sequential Adaptive Elastic Net Approach for Single-snapshot Source Localization	 
* __Short Title__: _Sequential Adaptive Elastic Net_ 
* __Authors__: Muhammad Naveed Tabassum and Esa Ollila
* __Maintainer__: Muhammad Naveed Tabassum
* __Date__: 19.03.2018

## Introduction 

The MATLAB package __*SAEN-LARS*__ provides an implementation of three algorithms proposed in our paper, titled: _"Sequential Adaptive Elastic Net Approach for Single-snapshot Source Localization"_.

In above paper and accordingly in this package, _sequential adaptive elastic net (SAEN)_ approach applies the complex-valued pathwise method in the weighted elastic-net framework, named as _c-PW-WEN_, sequentially by decreasing the sparsity level (order) from 3K to K in three stages. SAEN utilizes smartly chosen adaptive (i.e., data dependent) weights that are based on solutions obtained in the previous stage. c-PW-WEN algorithm computes the WEN solution paths for different values of EN tuning parameter and then selects the best solution. To achieve this in a computationally efficient way, we develop a homotopy method that is a complex-valued extension of the least angle regression and shrinkage (LARS) algorithm for weighted Lasso problem, which we refer to as _c-LARS-WLasso_. It is numerically cost effective and avoids an exhaustive grid-search over candidate values of the regularization parameter.

NOTE: The c-PW-WEN algorithm contains 
- [x] both Lasso and EN as special cases for unity weights. 
- [x] contains both adaptive Lasso and adaptive EN as special cases for data dependent weights.

## Demo | Example

The package contains a simple demo (Demo.mlx) that explains the usage of algorithms for direction-of-arrival (DoA) estimation with a uniform linear array (ULA) in compressed beamforming (CBF) application.
Moreover, an example (Example.m) for set-up 4 in the paper is also included in the package. To have repeatable results, the pseudorandom number generator settings, in terms of seed and type, are provided along with scenario data in the package as 'seed_data.mat'.

## Download | Usage

The _SAEN-LARS_ package contains following files:

__README__: This file.

__Functions__: Function files for implementation of three algorithms proposed in the paper.
> * _[saen.m](https://github.com/mntabassm/SAEN-LARS/blob/master/saen.m)_: The main function, sequential adaptive elastic net (SAEN) approach. 
> * _[cpwwen.m](https://github.com/mntabassm/SAEN-LARS/blob/master/cpwwen.m)_: Auxiliary function, called by the main function three times.
> * _[clarswlasso.m](https://github.com/mntabassm/SAEN-LARS/blob/master/clarswlasso.m)_: Auxiliary function for finding knots and respective solutions at found knots.

__Usage__: The files for the demo and an example.
> * _[Demo.mlx](https://github.com/mntabassm/SAEN-LARS/blob/master/Demo.mlx)_: A live script demo.
> * _[Example.m](https://github.com/mntabassm/SAEN-LARS/blob/master/Example.m)_: An example for DoA estimation with a ULA in CBF application.
> * _[seed_data.mat](https://github.com/mntabassm/SAEN-LARS/blob/master/seed_data.mat)_: Scenario data and pseudorandom number generator settings. 

Download the package and extract the files into a folder with “full control” permission.
Set the MATLAB home directory to the above folder. Thereafter, open either 'Demo.mlx' or 'Example.m' file in MATLAB and follow the steps.

[![DOI](https://zenodo.org/badge/125818623.svg)](https://zenodo.org/badge/latestdoi/125818623) 
