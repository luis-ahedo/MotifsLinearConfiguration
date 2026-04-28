function [LyaEx v]= GramSchmidtV1(LyaEx,v)

n = length(v);
E=eye(n,n);


for ii = 1:1:n
    % v(:,ii) = v(:,ii) / norm(v(:,ii));
    for jj = ii+1:1:n
        v(:,jj) = v(:,jj) - proj(v(:,ii),v(:,jj));
    end

    d=norm(v(:,ii));
    v(:,ii) = v(:,ii) / d;
    LyaEx=LyaEx+E(ii,:)*log(d);
end


% % %%% For debug
% dot(v(:,3),v(:,1))


    function w = proj(u,v)
        % This function projects vector v on vector u
        w = (dot(v,u) / dot(u,u)) * u;
    end

end