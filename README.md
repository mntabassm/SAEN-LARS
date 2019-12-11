# Software Package | _SAEN-LARS_    [![DOI](https://zenodo.org/badge/125818623.svg)](https://zenodo.org/badge/latestdoi/125818623)

* __Version__: 1.0.1 
* __GitHub Link__: https://github.com/mntabassm/SAEN-LARS

* __Title__: Sequential Adaptive Elastic Net Approach and Sparsity (or Model) Order Detection with an Application to Single-snapshot Source Localization	 
* __Short Title__: _Sequential Adaptive Elastic Net_ | _Sparsity (or Model) Order Detection_
* __Authors__: Muhammad Naveed Tabassum and Esa Ollila
* __Maintainer__: Muhammad Naveed Tabassum
* __Language__: MATLAB
* __Date__: 19.03.2018
* __Date (Last update)__: 20.06.2018

## Introduction 

The software package __*SAEN-LARS*__ provides an implementation (and examples) of algorithms proposed in our following papers: 

1. _"Sequential Adaptive Elastic Net Approach for Single-snapshot Source Localization"_ < https://doi.org/10.1121/1.5042363 >
2. _"Simultaneous Signal Subspace Rank and Model Selection with an Application to Single-snapshot Source Localization"_ < https://doi.org/10.23919/EUSIPCO.2018.8553171 >

In __*first paper*__ and accordingly in this package, _sequential adaptive elastic net (SAEN)_ approach applies the complex-valued pathwise method in the weighted elastic-net framework, named as _c-PW-WEN_, sequentially by decreasing the sparsity level (or order) from 3K to K in three stages. SAEN utilizes smartly chosen adaptive (i.e., data dependent) weights that are based on solutions obtained in the previous stage. The c-PW-WEN algorithm computes the WEN solution paths for different values of EN tuning parameter and then selects the best solution. To achieve this in a computationally efficient way, we develop a homotopy method that is a complex-valued extension of the least angle regression and shrinkage (LARS) algorithm for weighted Lasso problem, which we refer to as _c-LARS-WLasso_. It is numerically cost effective and avoids an exhaustive grid-search over candidate values of the regularization parameter.

NOTE: The c-PW-WEN algorithm contains both
- [x] Lasso and EN as special cases for unit weights. 
- [x] Adaptive Lasso and adaptive EN as special cases for data-dependent weights.

For __*second paper*__, we develop the _c-LARS-GIC_ method that is a two-stage procedure, where firstly precise values of the regularization parameter, called knots, at which a new predictor variable enters (or leaves) the active sets are computed in the Lasso solution path (using _c-LARS-WLasso_ with unit weights). Active sets provide a nested sequence of regression models and GIC then selects the best model using _c-LARS-GIC_. The sparsity order of the chosen model serves as an estimate of the model order.

## Demo | Example

The package contains a simple demo (Demo.mlx) that explains the usage of algorithms (of first paper) for direction-of-arrival (DoA) estimation with a uniform linear array (ULA) in compressed beamforming (CBF) application. Moreover, an example (Example.m) for setup 4 in the first paper is also included in the package. Another example (Example_GIC.m) in this package is for second simulation setup (i.e., Fig. 2) of the second paper, when the number of sensors in the ULA is n = 40.

NOTE: To have repeatable results, the pseudorandom number generator settings, in terms of seed and type, are provided along with scenarios data in the package as 'seed_data.mat' and 'seed_data_gic.mat'.

## Download | Usage

The _SAEN-LARS_ package contains following files:

__README__: This file.

__Functions__: Function files for implementation of algorithms proposed in both above mentioned papers.
> * _[saen.m](https://github.com/mntabassm/SAEN-LARS/blob/master/saen.m)_: The main function, sequential adaptive elastic net (SAEN) approach. 
> * _[cpwwen.m](https://github.com/mntabassm/SAEN-LARS/blob/master/cpwwen.m)_: Auxiliary function, called by the main function three times.
> * _[clarswlasso.m](https://github.com/mntabassm/SAEN-LARS/blob/master/clarswlasso.m)_: Auxiliary function for finding knots and respective solutions at found knots.
> * _[clarsgic.m](https://github.com/mntabassm/SAEN-LARS/blob/master/clarsgic.m)_: Auxiliary function for detecting the true sparsity (or model) order and estimating corresponding solution.

__Usage__: The files for the demo and examples.
> * _[Demo.mlx](https://github.com/mntabassm/SAEN-LARS/blob/master/Demo.mlx)_: A live script demo.
> * _[Example.m](https://github.com/mntabassm/SAEN-LARS/blob/master/Example.m)_: An example for DoA estimation with a ULA in CBF application.
> * _[seed_data.mat](https://github.com/mntabassm/SAEN-LARS/blob/master/seed_data.mat)_: Scenario data and pseudorandom number generator settings. 
> * _[Example_GIC.m](https://github.com/mntabassm/SAEN-LARS/blob/master/Example_GIC.m)_: An example for detection of sparsity (or model) order, i.e., the number of sources in CBF application.
> * _[seed_data_gic.mat](https://github.com/mntabassm/SAEN-LARS/blob/master/seed_data_gic.mat)_: Scenario data and pseudorandom number generator settings for application of c-LARS-GIC. 

Download the package and extract the files into a folder with “full control” permission.
Set the MATLAB home directory to the above folder. Thereafter, open either 'Demo.mlx', 'Example.m' or 'Example_GIC.m' file in MATLAB and follow the steps.

[![DOI](https://zenodo.org/badge/125818623.svg)](https://zenodo.org/badge/latestdoi/125818623) 
