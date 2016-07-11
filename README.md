# Matlab-WENO5
A matlab toolkit to calculate numerical differentiation using WENO5 scheme. Mainly for level set simulation.

## Install
Just Copy the code to your working directory or add them to your path.

``` addpath('X:\dev\GitHub\Matlab-WENO5');```

## Function usage

Currently the code only works with 2D. Extending to 3D is straight forward.

```function [ Dx,Dy,Dxf,Dxb,Dyf,Dyb ] = WENO5_2D( F,StepX,StepY )```

`F` is the 2D matrix you want to get WENO differentiation, `StepX,StepY` is grid size on `x` or `y` axis. Note that `x` here corresponds to the second subscript in a matlab matrix.

```
%Caculate numerical differentiation using WENO5 scheme, and choose one from
%backwrd/forward differentiation using Godunov's scheme. In this file
%"x-axis" is left-right, i.e. the second dimmendion in a matlab matrix, and
%"y-axis" is in top-down direction
%   F: input scalr field
%   StepX, StepY: Grid size on the two directions
%   Dx: output differentiation value on x-axis, according to Godunov's scheme
%   Dy: output differentiation value on y-axis, according to Godunov's scheme
%   Dxf,Dxb,Dyf,Dyb: WENO5 backwrd/forward differentiation on x/y axis
%       "f" maens forward and "b" maens backward
```

## Further work

Extend the code to 3D, add test or example code. Any contribution will be welcomed.