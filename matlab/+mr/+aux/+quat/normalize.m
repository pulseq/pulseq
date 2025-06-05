function q = normalize(q)
%NORMALIZE normalizes quaternion q (or array of quaternions)

n2=sum(q'.^2)';
if any(n2>0)
    if length(n2)>1
        ninv=n2(n2>0).^-0.5;
        q(n2>0,:)=q(n2>0,:).*ninv(:,ones(1,4));
    else
        ninv=n2^-0.5;
        q=q*ninv;
    end
end
