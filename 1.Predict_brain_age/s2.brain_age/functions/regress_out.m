
function [r, b, stats] = regress_out(Y, conf)
% Y =[n,1], conf=[n,q], there are q confound variables
    [M,N] = size(Y);
    Y = reshape(Y, [length(Y),1]);
    X = [ones(size(conf,1),1), conf];
    [b,bint,r,rint,stats] = regress(Y, X);
    r = reshape(r, [M,N]);
end
