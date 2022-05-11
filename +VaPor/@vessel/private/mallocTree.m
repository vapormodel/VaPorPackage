function [treeOut] = mallocTree(treeIn, iter)
%MALLOCTREE Allocates memory for the new tree.
%   This function creates the new matrix to hold the expanded tree.
%   Author:         Luke Fulford - luke.fulford@ed.ac.uk
%   Last Modified:  29/01/2020

numVessel = size(treeIn,1); % Number of points comprising the current vessel tree.
treeOut = zeros(numVessel+2*iter,7); % Preallocate maximum possible size of new vessel tree (current size + 2*no_iter).
treeOut(1:numVessel,:) = treeIn; % Insert original tree into new tree to start.
end

