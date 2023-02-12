function g = dom_g(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn)

n = length(x);
g = zeros(n,1);
ksum = zeros(n,1);
for k=1:n
    if(cnn(invindex(k)) == 0)
        ksum(k) = -alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1);
    end
end


for i=1:n
    g(i) = calc_obj*(1+delta) - calc_obj*2*delta*x(i) - beta(i)/x(i) + gamma(i)/(1-x(i));
    %disp(['i = ' num2str(i)]);
    for k=index(cn(i,1:degrees(i)+1))'
        if(k == -1)
            continue;
        end
        if(cnn(invindex(k)) == 0)
    %    k
            %g(i) = g(i) - alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1);
            g(i) = g(i) + ksum(k);
        end  
    end
    %pause
end