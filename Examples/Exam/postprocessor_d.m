function [z, M] = postprocessor_d(state)
R0 = 286.9;
M0 = 29.97;
% compute compressibilty
z = state.p / (state.r*R0*state.T);

% compute molecular weight
M = M0/z;
