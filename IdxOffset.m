function [ r ] = IdxOffset( crd, crdLim, offset )
%Calculate bounded coordinate offset. make sure matrix subset never go out
%of matrix size
%   crd: coordinate to be offset
%   crmLim: upper bound of coordinate (low bound is always 1)

if (crd + offset < 1)
    r = 1;
    return
end
if (crd + offset > crdLim)
    r = crdLim;
    return;
end

r = crd  + offset;

end