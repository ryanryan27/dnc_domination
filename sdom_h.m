function h = sdom_h(x, nbrs, cn, degrees, two, two_degree, ind2ij,ij2ind,alpha, beta,gamma)

n = size(two_degree,1);

vars = length(x);

h = zeros(vars,vars);


for i=1:n
    
    %d2fdx2 rows 3,4,5,9
    h(i,i) = h(i,i) + alpha(1)/(x(i)^2) + alpha(2)/((1-x(i))^2) - 2*beta(1) + 2*beta(2);
    
    for q=1:n
       for k = intersect(cn(i,1:degrees(i)+1),cn(q,1:degrees(q)+1))
          
           for j = two(k,1:two_degree(k))
                h(i,q) = h(i,q) + gamma(3)/((sum(x(cn(k,1:degrees(k)+1))) - sum(y(intersect(nbrs(k,1:degrees(k)),nbrs(j,1:degrees(j))),j)) - 1)^2);
               
           end
           
           
           
           
       end
        
    end
    
    
    
    
    
    
    
    for k = cn(i,1:degrees(i)+1)
%         inner = 0;
%         for m = two(k,1:two_degree(k))
%             %inner calc for d2fdxdx row 8
%             inner = sum(x(cn(k,1:degrees(k)+1)));
%             inner = inner - sum(y(intersect(nbrs(k,1:degrees(k)),nbrs(m,1:degrees(m))),m));
%             
%         end
        
        for j = cn(k,1:degrees(k)+1)
            %d2fdxdx d2fdx2 row 7
            h(i,j) = h(i,j) + alpha(3)/((sum(x(cn(k,1:degrees(k)+1)))-1)^2);
%             
%             %d2fdxdx d2fdx2 row 8
%             h(i,j) = h(i,j) + gamma(3)/((inner -1)^2);
%             
        end
        
        
    end
    
    for j = nbrs(i,1:degrees(i))
        
        %d2fdx2 row 6
        h(i,i) = h(i,i) + gamma(2)/((x(i)-y(i,j))^2);
        
        %d2fdy2 row 2
        h(ij2ind(i,j), ij2ind(i,j)) = h(ij2ind(i,j), ij2ind(i,j)) + gamma(1)/(y(i,j)^2);
        
        %d2fdydx row 6
        h(ij2ind(i,j), i) = h(ij2ind(i,j), i) - gamma(2)/((x(i)-y(i,j))^2);
        h(i, ij2ind(i,j)) = h(i, ij2ind(i,j)) - gamma(2)/((x(i)-y(i,j))^2);
        
        
        %d2fdy2 row 6
        h(ij2ind(i,j), ij2ind(i,j)) = h(ij2ind(i,j), ij2ind(i,j)) + gamma(2)/((x(i)-y(i,j))^2);
        
        %
        %d2fdydx row 9
        h(ij2ind(i,j), j) = h(ij2ind(i,j), j) + 2*beta(2);
        h(j,ij2ind(i,j)) = h(j,ij2ind(i,j)) + 2*beta(2);
        
  
        %d2fdy2 row 9
%         h(ij2ind(i,j), ij2ind(i,j)) = h(ij2ind(i,j), ij2ind(i,j)) + 2*beta(2);
        
        for k = nbrs(j,1:degrees(j))
            
        
            %d2fdydy row 9
            h(ij2ind(i,j), ij2ind(k,j)) = h(ij2ind(i,j), ij2ind(k,j)) + 2*beta(2);
            
        end
        
        for k = two(j,1:two_degree(j))
            
            if(~ismember(k,nbrs(i,1:degrees(i))))
                continue;
            end
            
            inner = sum(x(cn(k,1:degrees(k)+1)));
            inner = inner - sum(y(intersect(nbrs(k,1:degrees(k)),nbrs(j,1:degrees(j))),j));
            
            %if(ismember(i,nbrs(k,1:degrees(k))))
                %d2fdy2 row 8
            h(ij2ind(i,j), ij2ind(i,j)) = h(ij2ind(i,j), ij2ind(i,j)) + gamma(3)/((inner -1)^2);
            %end
            
            
            
            for l = cn(k,1:degrees(k)+1)
                %d2fdydx row 8
                h(ij2ind(i,j), l) = h(ij2ind(i,j), l) - gamma(3)/((inner-1)^2);
                h(l,ij2ind(i,j)) = h(l,ij2ind(i,j)) - gamma(3)/((inner-1)^2);
                
            end
            
            for p = intersect(nbrs(k,1:degrees(k)),nbrs(j,1:degrees(j)))
                if(p~=i)
                    %d2fdydy row 8
                    h(ij2ind(i,j), ij2ind(p,j)) = h(ij2ind(i,j), ij2ind(p,j)) + gamma(3)/((inner-1)^2);
                end
            end
            
            
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