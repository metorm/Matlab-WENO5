function [ D ] = CombineLR( derivL,derivR,S )
%Gudunov's method to choose one from backward/forward differentiation
%   derivL: matrix of backward differentiation
%   derivR: matrix of forward differentiation
%   S: source scalar field
%   D: return value

if(sum(abs(size(S)-size(derivL))) + sum(abs(size(S)-size(derivL))) > 0)
    error('Size of input arrays are different!')
end

flowL = ((S .* derivR <= 0) & (S .* derivL <= 0));

% Both directions agree that flow is to the right.
flowR = ((S .* derivR >= 0) & (S .* derivL >= 0));

% Diverging flow; entropy condition requires choosing deriv = 0
%   (so we don't actually have to calculate this term).
%flow0 = ((S .* derivR >  0) & (S .* derivL <  0));

% Converging flow, need to check which direction arrives first.
flows = ((S .* derivR <  0) & (S .* derivL >  0));
if(any(flows(:)))
  conv = find(flows);
  s = zeros(size(flows));
  s(conv) = S(conv) .* (abs(derivR(conv)) - abs(derivL(conv))) ...
            ./ (derivR(conv) - derivL(conv));

  % If s == 0, both directions arrive at the same time.
  %   Assuming continuity, both will produce same result, so pick one.
  flowL(conv) = flowL(conv) | (s(conv) < 0);
  flowR(conv) = flowR(conv) | (s(conv) >= 0);
end

D = derivL .* flowR + derivR .* flowL;

if(sum(sum(isnan(D)))>0)
    % if this error is triggered, set break point here to debug
    error('NaN found in array. Something wrong.');
end

end