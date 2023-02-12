function g = sdom_g(x, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma)


 n = size(two_degree,1);

 vars = length(x);
 
g = zeros(vars,1);

for i=1:n
    
    %dfdx rows 1,3,4,5
    g(i) = 1 - alpha(1)/x(i) + alpha(2)/(1-x(i)) + beta(1)*(1-2*x(i));
    
    %dfdx row 9
    g(i) = g(i) + 2*beta(2)*(x(i) + sum(y(nbrs(i,1:degrees(i)),i))-1);
    
    
    for h=cn(i,1:degrees(i)+1)
        %dfdx row 7
        g(i) = g(i) - alpha(3)/(sum(x(cn(h,1:degrees(h)+1)))-1);
        
        %dfdx row 8
        for j=two(h,1:two_degree(h))
            
            inner = sum(x(cn(h,1:degrees(h)+1)));       
            inner = inner - sum(y(intersect(nbrs(j,1:degrees(j)),nbrs(h,1:degrees(h))),j));
            g(i) = g(i) - gamma(3)/(inner - 1);
            
        end
    end
    
    
    for j = nbrs(i,1:degrees(i))
        
        %dfdx row 6
        g(i) = g(i) - gamma(2)/(x(i) - y(i,j));
        
        %dfdy rows 2,6
        g(ij2ind(i,j)) = gamma(2)/(x(i) - y(i,j))-1*gamma(1)/y(i,j);  
        
        
        %dfdy row 9 
        g(ij2ind(i,j)) = g(ij2ind(i,j)) + 2*beta(2)*(x(j)+sum(y(nbrs(j,1:degrees(j)),j))-1);
        
        %dfdy row 8
        for k=two(j,1:two_degree(j))
            
            if(~ismember(i,nbrs(k,1:degrees(k)))) 
                continue;
            end
            
            inner = sum(x(cn(k,1:degrees(k)+1)));       
            inner = inner - sum(y(intersect(nbrs(k,1:degrees(k)),nbrs(j,1:degrees(j))),j));
            
            g(ij2ind(i,j)) = g(ij2ind(i,j)) + gamma(3)/(inner - 1);
            
        end
        
        
    end
    
end


    function val = y(a,b)
        
        if(length(a) > length(b))
            b = b*ones(1,length(a));
        end
        if(length(b) > length(a))
            a = a*ones(1,length(b));
        end

        val = sum(x(diag(ij2ind(a,b))));
        
    end



end