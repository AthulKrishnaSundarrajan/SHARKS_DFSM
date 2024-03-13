function [A,B,C,D] = LTI_function(x,Ts,ns,nc)


% reshape
x = reshape(x,[ns/2,ns+nc]);

% extract right elements
B_par = x(:,1:nc);
A_par = x(:,nc+1:nc+ns);

A = [zeros(ns/2),eye(ns/2);A_par];


B = [zeros(ns/2,nc);B_par];
C = A_par; 
D = B_par; 

end