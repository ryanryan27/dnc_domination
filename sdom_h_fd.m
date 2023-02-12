function H = sdom_h_fd(x,nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma)


h=0.0001;
    
vars = length(x);


H = zeros(vars,vars);

for i = 1:vars
    
    for j = 1:vars
    
        if(i == j)
           
            H(i,j) = (-1*f(x+2*h*e_vector(i,vars)) + 16*f(x+h*e_vector(i,vars)) - 30*f(x) + 16*f(x-h*e_vector(i,vars)) - f(x - 2*h*e_vector(i,vars)))/(12*h^2);
            
        else
            
            H(i,j) = (f(x+h*e_vector(i,vars) + h*e_vector(j,vars)) - f(x+h*e_vector(i,vars) - h*e_vector(j,vars)) - f(x-h*e_vector(i,vars) + h*e_vector(j,vars)) + f(x - h*e_vector(i,vars) - h*e_vector(j,vars)))/(4*h^2);
            
            
        end
        
        
        
    end
    
end




    function val = f(in)
        val = sdom_f(in, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind, alpha, beta,gamma);
        
    end


end