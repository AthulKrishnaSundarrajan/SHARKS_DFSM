function [A,B] = mufunc_AB(par,ns,nc)

% reshape
par = reshape(par,[ns/2,ns+nc]);

% extract right elements
B_par = par(:,1:nc);
A_par = par(:,nc+1:ns);

A = [zeros(ns/2),eye(ns/2);A_par];
B = [zeros(nc);B_par];




end