function [ Dx,Dy ] = BackwardFirstOrderDiff( F,StepX,StepY )
%Calculate backward first order differentiation for filed F. dx_n = (x_n - x_(n-1))/dx
%   F: input scalr field
%   StepX, StepY: Grid size on the two directions
%   Dx: output differentiation value on x-axis (left-right, second dimmendion)
%   Dy: output differentiation value on y-axis (top-down, first dimension)

if(StepX==0 || StepY==0 || size(F,1)<2 || size(F,2)<2 )
    error('Invalid input!');
end

Dy = zeros(size(F));
Dx = zeros(size(F));

Dx(:,2:end) = diff(F,1,2)/StepX;
Dx(:,1) = Dx(:,2);

Dy(2:end,:) = diff(F)/StepY;
Dy(1,:) = Dy(2,:);

if(sum(sum(isnan(Dx))) + sum(sum(isnan(Dy))) >0)
    % if this error is triggered, set break point here to debug
    error('NaN found in array. Something wrong.');
end

end