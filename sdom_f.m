function f = sdom_f(x, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma)
   

    n = size(two_degree,1);
    
    %row 1
    f = sum(x(1:n));
    
    for i=1:n
        %row 3
        f = f - alpha(1)*log(x(i));
        
        %row 4
        f = f - alpha(2)*log(1-x(i));
        
        %row 5
        f = f + beta(1)*(x(i)-x(i)^2);
        
        %row 7
        f = f - alpha(3)*log(sum(x(cn(i,1:degrees(i)+1)))-1);
        
        %row 9
        f = f + beta(2)*(x(i) + sum(y(nbrs(i,1:degrees(i)),i))-1)^2;
        
        for j=nbrs(i,1:degrees(i))
            %row 2
            f = f - gamma(1)*log(y(i,j));

            %row 6
            f = f - gamma(2)*log(x(i)-y(i,j));
            if(gamma(2)*log(x(i)-y(i,j)) < -10)
                disp(['contribution: ' gamma(2)*log(x(i)-y(i,j))])
                pause
            end
            
        end
        
        
        for j=two(i,1:two_degree(i))
            
            %row 8
            inner = sum(x(cn(j,1:degrees(j)+1)));       
            inner = inner - sum(y(intersect(nbrs(i,1:degrees(i)),nbrs(j,1:degrees(j))),i));
            f = f - gamma(3)*log(inner - 1);
            
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