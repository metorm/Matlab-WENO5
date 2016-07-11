function [ Dx,Dy,Dxf,Dxb,Dyf,Dyb ] = WENO5_2D( F,StepX,StepY )
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

% calculate first order backward differentiation. xf_n = xb_(n+1)
[BackwardDx,BackwardDy] = BackwardFirstOrderDiff( F,StepX,StepY );%#ok<*PFBNS>

if(sum(sum(isnan(BackwardDx))) + sum(sum(isnan(BackwardDy))) >0)
    % if this error is triggered, set break point here to debug
    error('NaN found in array. Something wrong.');
end

[LimY,LimX] = size(F);

% pre-allocate
Dxf=zeros(size(F));Dxb=Dxf;Dyf=Dxf;Dyb=Dxf;

% for-each calculate
for xi = 1:LimX
    parfor yi = 1:LimY
        
        % frequently used value
        DiffBackwardValueX = BackwardDx(yi,xi);
        DiffBackwardValueY = BackwardDy(yi,xi);
        
        Dxb(yi,xi) = WENO5Kernel(...
            BackwardDx(yi,IdxOffset(xi,LimX,-2)),...
            BackwardDx(yi,IdxOffset(xi,LimX,-1)),...
            DiffBackwardValueX,...
            BackwardDx(yi,IdxOffset(xi,LimX,+1)),...
            BackwardDx(yi,IdxOffset(xi,LimX,+2)));
        
        Dxf(yi,xi) = WENO5Kernel(...
            BackwardDx(yi,IdxOffset(xi,LimX,+3)),...
            BackwardDx(yi,IdxOffset(xi,LimX,+2)),...
            BackwardDx(yi,IdxOffset(xi,LimX,+1)),...
            DiffBackwardValueX,...
            BackwardDx(yi,IdxOffset(xi,LimX,-1)));
        
        Dyb(yi,xi) = WENO5Kernel(...
            BackwardDy(IdxOffset(yi,LimY,-2),xi),...
            BackwardDy(IdxOffset(yi,LimY,-1),xi),...
            DiffBackwardValueY,...
            BackwardDy(IdxOffset(yi,LimY,+1),xi),...
            BackwardDy(IdxOffset(yi,LimY,+2),xi));
        
        Dyf(yi,xi) = WENO5Kernel(...
            BackwardDy(IdxOffset(yi,LimY,+3),xi),...
            BackwardDy(IdxOffset(yi,LimY,+2),xi),...
            BackwardDy(IdxOffset(yi,LimY,+1),xi),...
            DiffBackwardValueY,...            
            BackwardDy(IdxOffset(yi,LimY,-1),xi));        
    end
end

if(sum(sum(isnan(Dxf))) + sum(sum(isnan(Dxb))) + sum(sum(isnan(Dyb))) ...
    + sum(sum(isnan(Dyf))) > 0)
    % if this error is triggered, set break point here to debug
    error('NaN found in array. Something wrong.');
end

% combine backward and forward differentiation
Dx=CombineLR(Dxb,Dxf,F);
Dy=CombineLR(Dyb,Dyf,F);
end