function mse = calculate_mse(X1,X2)

    % find the shape
    [nt,nX] = size(X1);
    
    % initialize storage array
    mse = zeros(nX,1);
    
    % loop through the second dimesnion and calculate MSE
    for ix = 1:nX

        mse(ix) = sum((X1(:,ix)-X2(:,ix)).^2)/nt;
        
    end


end