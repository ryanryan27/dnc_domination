function [x,t,k] = domopt(A,seed)
tic;
%rand('state',seed);
rng(seed);
displays = 0;

fixcutoff = 0.99;

n = size(A,1);
orig_n = n;
calc_obj = 0;
alpha = 1*ones(n,1);
beta = 0.1*ones(n,1);
gamma = 0.1*ones(n,1);
delta = 0.1;



%x = 0.5*ones(n,1);
x = 0.5+0.5*rand(n,1);
degrees = sum(A);
cn = zeros(n,max(degrees)+1);
for i=1:n
    count = 0;
    for j=1:n
        if(A(i,j) || i == j)
            count = count + 1;
            cn(i,count) = j;
        end
    end
end

ind = 0;
newindex = 1;
ijintind(newindex) = 1;
for i=1:n
    for j=i:n
        count = 0;
        for k=cn(i,1:degrees(i)+1)
            if((A(j,k) || j == k))
                count = count + 1;
                ind = ind + 1;
                ijint(ind) = k;                
            end
        end
        if(count > 0)
            ijinti(newindex) = i;
            ijintj(newindex) = j;
            newindex = newindex + 1;
            ijintind(newindex) = length(ijint)+1;            
        end
    end
end



rows = sqrt(n);
cols = sqrt(n);
coords = zeros(n,2);

count = 0;
for i=1:rows
    for j=1:cols
        count = count + 1;
        coords(count,[1 2]) = [j,rows+1-i];
    end
end




index = [1:n]';
invindex = [1:n]';
cnn = zeros(n,1);


f = dom_f(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
g = dom_g(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
H = dom_h(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn,A,ijint,ijintind,ijinti,ijintj);

DDcount = 0;
DNCcount = 0;
dtype = '';
iterations = 0;
while 1
    [a,b] = eig(H);
    dnc = a(:,1);
    se = b(1,1);
    stopiterating = 0;
    while(norm(g) > 1e-4 || se < -1e-4)
        xreal = 2*ones(orig_n,1);
%         invindex
%         x
        xreal(invindex) = x;
        if(displays >= 2)
            domvis(A,coords,xreal);
            pause(0.01);
        end
        if(strcmp(dtype,'DNC'))
            
            %pause;
            
        end
        %pause(0.1)
        if(se < -1e-4)
            %disp('DNC');
            DNCcount = DNCcount + 1;
            d = dnc;
            dd = -inv(H)*g;
            dtype = 'DNC';
            
            %d'
            %pause
        else
            %disp('DD');
            DDcount = DDcount + 1;
            d = zeros(size(H,1),1);
            dd = -inv(H)*g;
            dtype='DD';
        end

        if(g'*d > 1e-4)
            d = -d;
        end
        if(g'*dd > 1e-4)
            dd = -dd;
        end
        
        d = d + dd;
        
        step = -1;
        for i=1:n
            if(d(i) > 0)
                if(step >= 0.9*(1-x(i))/d(i) || step == -1)
                    step = 0.9*(1-x(i))/d(i);
                end                   
            else
                if(step >= 0.9*(-x(i))/d(i) || step == -1)
                    step = 0.9*(-x(i))/d(i);
                end                   
            end
        end

        for i=1:n
            if(cnn(invindex(i)) == 0)
                dsum = sum(d(index(cn(i,1:degrees(i)+1))));
                if(dsum < 0)
                    if(step >= 0.9*(1 - sum(x(index(cn(i,1:degrees(i)+1)))))/dsum)
                        step = 0.9*(1 - sum(x(index(cn(i,1:degrees(i)+1)))))/dsum;
                    end
                end
            end
        end
        iters = 0;
        while(dom_f(x+step*d,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn) >= f)
            step = step/2;
            
            iters = iters + 1;
            
            if(iters > 1000)
                x = ones(orig_n,1);
                t = toc;
                k = orig_n;
                disp(['got stuck with seed: ' num2str(seed)]);
                return  
            end
            
        end
           
        
        
        
        
        %step
        x = x + step*d;
        
        %%% CHECK FOR FIXING TO 1
        
        checkfix = 1;
        while(checkfix && calc_obj)
            checkfix = 0;
            for i=1:n
                if(x(i) > fixcutoff)
                    if(displays >=1)
                        disp(['FIXING VARIABLE ' num2str(i) ' (originally var ' num2str(invindex(i)) ')']);
                    end
                    checkfix = 1;
                    index(invindex(i)) = -1;
                    for j=invindex(i)+1:orig_n
                        if(index(j) ~= -1)
                            index(j) = index(j) - 1;
                        end
                    end
                    invindex = find(index ~= -1);
%                     invindex = zeros(n-1,1);
%                     invcount = 0;
%                     for j=1:length(index)
%                         invcount = invcount + 1;
%                         if(index(j) ~= -1)
%                             invindex(index(j)) = j;
%                         end
%                     end
                    cnn(cn(i,1:degrees(i)+1)) = 1;
                    x = x([1:i-1 i+1:n]);
                    alpha = alpha([1:i-1 i+1:n]);
                    beta = beta([1:i-1 i+1:n]);
                    gamma = gamma([1:i-1 i+1:n]);
                    cn = cn([1:i-1 i+1:n],:);
                    degrees = degrees([1:i-1 i+1:n]);
                    n = n - 1;
                    break;
                end
            end
        end
        
        
        
        f = dom_f(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
        
        iterations = iterations + 1;
        realf(iterations) = sum(x) + (orig_n - length(x));
        [~,ind_far] = min(abs(x-0.5));
        if(displays)
            displayer(iterations, dtype, 'f',f,2, 'real f',realf(iterations),2, 'se',se,2,  'norm(g)',norm(g),2, 'g''*d',g'*d, 2, 'step', step, 2, 'max_x', max(x), 2, 'min_x', min(x),2, 'farthest_x', x(ind_far),2);
        end
        
        g = dom_g(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
        H = dom_h(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn,A,ijint,ijintind,ijinti,ijintj);
        %H
        f;
        [a,b] = eig(H);
        dnc = a(:,1);
        se = b(1,1);
        stopiterating = 1;
        for i=1:n
            if(x(i) < 0.99 && x(i) > 0.01)
                stopiterating = 0;
                break;
            end
        end
        if(stopiterating)
            break;
        end
    end
    if(stopiterating)
        x = round(x);
        if(displays >=1)
            disp(['FINAL NUMBER OF GUARDS : ' num2str(sum(x) + sum(index == -1))]);
        end
        break;
    end
    alpha = 0.5*alpha;
    beta = 0.5*beta;
    gamma = 0.5*gamma;
    delta = 1.5*delta;
    calc_obj = 1;
    if(displays >= 1)
        disp('REDUCING PARAMETERS')
        param_reduce(iterations) = 1;
    end
    f = dom_f(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
    g = dom_g(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn);
    H = dom_h(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn,A,ijint,ijintind,ijinti,ijintj);
    
    
end
% DNCcount
% DDcount
x = (round(xreal) > 0);
k = sum(x);
t=toc;
if(displays >= 1)
    figure(2)
    plot(1:iterations,realf)
    hold on;
    plot(find(param_reduce),realf(find(param_reduce)),'x')
    hold off;
end