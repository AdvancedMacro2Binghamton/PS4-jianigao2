%%%%%%%%%parameters%%%%%%%%
alpha=1/3;
beta=0.99;
s=2;
delta=0.025;
k=30;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=5;
rho=0.5;
sigma=0.2;
d=3;

[z,zprob] =TAUCHEN(m,rho,sigma,d);
z=exp(z);
pi=zprob^1000;
pz=pi(1,:);
ns=pz*z;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_min=0;
a_max=300;
a_num=100;
a=linspace(a_min,a_max,a_num);

r=alpha*k_guess^(alpha-1)*nd^(1-alpha)+1-delta;
w=(1-alpha)*k_guess^alpha*nd^(-alpha);

%consumption fn
cons=bsxfun(@minus,r*a',a);
cons=bsxfun(@plus,cons,permute(w*z,[2 3 1]));
%return fn
ret=cons.^(1-s)/(1-s);
ret(cons<0)=-Inf;
%value fn and policy fn
v_guess=zeros(m,a_num);

        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
v_tol=1;

while v_tol>1e-8
    
vf=bsxfun(@plus,ret,permute(beta*zprob*v_guess,[3,2,1]));
% CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
[vfn,pol_indx]=max(vf,[],2);
pol_indx=shiftdim(pol_indx,2);  
pol_fn=a(pol_indx);  
vfn=shiftdim(vfn,2);

Q = makeQmatrix(pol_indx,zprob);
u=bsxfun(@minus,r*a,pol_fn);
u=u+repmat(z*w,[1 a_num]);
ret=u.^(1-s)/(1-s);
ret(u<0)=-Inf;

u=ret(:);
w_vec=v_guess(:);
 for ii=1:k
    w_new=u+beta*Q*w_vec;
    w_vec=w_new;
 end

v_guess=w;

v_tol=max(abs(vfn(:)-v_guess(:)));
end

%%%%%%%%%question 6
plot(a,pol_fn(1,:),a,pol_fn(2,:),a,pol_fn(3,:),a,pol_fn(4,:),a,pol_fn(5,:)),legend('z1','z2','z3','z4','z5');

%gini
p=[Mu(1,:);Mu(2,:);Mu(3,:);Mu(4,:);Mu(5,:)];
wealth=[a;a;a;a;a];
wg=gini(p,wealth,true);
title(['wealth gini index=',num2str(wg)]);

