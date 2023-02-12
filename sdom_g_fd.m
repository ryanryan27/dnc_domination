function g = sdom_g_fd(x, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma)


h = 0.0001;
n = length(two_degree);
vars = length(x);
g = zeros(vars,1);

for i = 1:vars
    g(i) = (f(x+ h*e_vector(i,vars)) - f(x-h*e_vector(i,vars)))/(2*h);
    
    
end


    function val = f(in)
        val = sdom_f(in, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma);
        
    end

end