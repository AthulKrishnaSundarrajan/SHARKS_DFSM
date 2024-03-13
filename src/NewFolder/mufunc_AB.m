function O = mufunc_AB(x,ns,nc,inputs,state_dx,R,gamma)

% reshape
x_ = x;
x = reshape(x,[ns/2,ns+nc]);

% extract right elements
B_par = x(:,1:nc);
A_par = x(:,nc+1:nc+ns);

A = [zeros(ns/2),eye(ns/2);A_par];
B = [zeros(ns/2,nc);B_par];


LM = [A,B];

LM = LM';

error = norm((state_dx - inputs*LM).^2);

regularization = gamma*(x_'*R*x_);



O = error + regularization;


end