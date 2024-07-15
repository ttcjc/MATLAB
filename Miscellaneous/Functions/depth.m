function D = depth(X)
%DEPTH Number of pages in an array.
%   D = DEPTH(X) returns the number of pages in the array X.  DEPTH(X) is
%   equivalent to SIZE(X,3).
%
%   See also HEIGHT, WIDTH, SIZE, NUMEL.

D = size(X,3);