function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1)
% function ulimit = SlopeLimitLin(ul,xl,vm1,v0,vp1);
% Purpose: Apply slopelimited on linear function ul(Np,1) on x(Np,1)
%          (vm1,v0,vp1) are cell averages left, center, and right

Variables;

% Compute various geometric measures
ulimit = ul; h = xl(P+1,:)-xl(1,:); 
x0 = ones(P+1,1)*(xl(1,:) + h/2);

hN = ones(P+1,1)*h;

% Limit function
ux = (2./hN).*(Dr*ul);

ulimit = ones(P+1,1)*v0+(xl-x0).*(ones(P+1,1)*minmod([ux(1,:);(vp1-v0)./h;(v0-vm1)./h]));
return
