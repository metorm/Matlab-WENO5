function [ D ] = CombineLR( L,R,F )
%Gudunov's method to choose one from backward/forward differentiation
%   L: matrix of backward differentiation
%   R: matrix of forward differentiation
%   F: source scalar field
%   D: return value

if(sum(abs(size(F)-size(L))) + sum(abs(size(F)-size(R))) > 0)
    error('Size of input arrays are different!')
end

S = (F>0)*2-1;

SR = S.*R;
SL = S.*L;

% For Gudunov's method, check characteristic directions according to left and right derivative approximations.
% Both directions agree that flow is to the left.
flowL = ((SR <= 0) & (SL <= 0));

% Both directions agree that flow is to the right.
flowR = ((SR >= 0) & (SL >= 0));

%Converging flow, need to check which direction arrives first.
flows = ((SR < 0) & (SL >  0));

for i = 1:numel(flows)
    if (flows(i))
        s = S(i) * (abs(R(i)) - abs(L(i))) / (R(i) - L(i));
        % If s == 0, both directions arrive at the same time. Assuming continuity, both will produce same result, so pick one.
        flowL(i) = flowL(i) | (s < 0);
        flowR(i) = flowR(i) | (s >= 0);
    end
end

D = L.*flowR + R.*flowL;

end