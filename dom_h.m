function H = dom_h(x,cn,degrees,calc_obj,alpha,beta,gamma,delta,index,invindex,cnn,A,ijint,ijintind,ijinti,ijintj)

n = length(x);
% H = zeros(n);
% 
% for i=1:n
%     kindices = index(cn(i,1:degrees(i)+1))';
%     for j=1:n
%         if(j == i)
%             H(i,i) = -2*delta*calc_obj + beta(i)/x(i)^2 + gamma(i)/(1-x(i))^2;
%         end
%         for k=kindices
%             if(k == -1)
%                 continue;
%             end
%             if(cnn(invindex(k)) == 0)
%                 if(j == i)
%                     H(i,i) = H(i,i) + alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1)^2;
%                 else
%                     for l=1:degrees(j)+1
%                         if(index(cn(j,l)) == k)
%                             H(i,j) = H(i,j) + alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1)^2;
%                         end
%                     end
%                 end
%             end
%         end
%     end
% end

H = zeros(n);

for i=1:n
    %kindices = index(cn(i,1:degrees(i)+1))';
    H(i,i) = -2*delta*calc_obj + beta(i)/x(i)^2 + gamma(i)/(1-x(i))^2;
end

ksum = zeros(n,1);
for k=1:n
    if(cnn(invindex(k)) == 0)
        ksum(k) = alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1)^2;
    end
end

for place=1:length(ijintind)-1
    i = index(ijinti(place));
    j = index(ijintj(place));
    if(i == -1 || j == -1)
        continue;
    end
    for ind = ijintind(place):ijintind(place+1)-1
        k = index(ijint(ind));
        if(k == -1)
            continue;
        end
        if(cnn(invindex(k)) == 0)
            %H(i,j) = H(i,j) + alpha(k)/(sum(x(index(cn(k,1:degrees(k)+1))))-1)^2;
            H(i,j) = H(i,j) + ksum(k);
            H(j,i) = H(i,j);
        end
    end
end


% if(abs(sum(sum(abs(H - H)))) > 1e-4)
%     abs(sum(sum(abs(H - H))))
%     H
%     H
%     pause
% end
