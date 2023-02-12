function f = dom_f(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn)

n = length(x);
f = calc_obj*(1+delta)*sum(x);
for i=1:n
    f = f - calc_obj*delta*x(i)^2 - beta(i)*log(x(i)) - gamma(i)*log(1-x(i));
    if(cnn(invindex(i)) == 0)
        f = f - alpha(i)*log(sum(x(index(cn(i,1:degrees(i)+1))))-1);
    end
end