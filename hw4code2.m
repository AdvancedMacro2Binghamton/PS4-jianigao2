%%%%%%%%%parameters%%%%%%%%
alpha=1/3;
beta=0.99;
s=2;
delta=0.025;


%%%%%question 3%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m=5;
rho=0.5;
sigma=0.2;
d=3;

[z,zprob] =TAUCHEN(m,rho,sigma,d);
z=exp(z);
pi=zprob^1000;
pz=pi(1,:);
ns=pz*z;
%%%%%%%question 4%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_min=0;
a_max=5;
a_num=10;
a=linspace(a_min,a_max,a_num);
%%%%%%question 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nd=ns;
k_guess=2;


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
      
        vf=bsxfun(@plus,ret,permute(beta*zprob*v_guess,[3,2,1]));
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn,pol_indx]=max(vf,[],2);
        
    pol_fn=a(pol_indx); 
    a=interp(pol_fn(:,:,1),10);
    b=interp(pol_fn(:,:,2),10);
    c=interp(pol_fn(:,:,3),10);
    d=interp(pol_fn(:,:,4),10);
    e=interp(pol_fn(:,:,5),10);
    pol_indx=permute(pol_indx, [3 1 2]);
   
Q = makeQmatrix(pol_indx,zprob)
    
cons_new=bsxfun(@minus,r*a',permute(a(pol_indx),[2 3 1]));
cons_new=permute(cons_new,[3 1 2]);
cons_new=bsxfun(@plus,cons_new,w*z);   
ret_new=cons_new.^(1-s)/(1-s);
ret_new(cons_new<0)=-Inf;
r=[ret_new(:,1);ret_new(:,2);ret_new(:,3);ret_new(:,4);ret_new(:,5);ret_new(:,6);ret_new(:,7);ret_new(:,8);ret_new(:,9);ret_new(:,10)];
v_n=permute(vfn,[3 1 2]);
v=[v_n(:,1);v_n(:,2);v_n(:,3);v_n(:,4);v_n(:,5);v_n(:,6);v_n(:,7);v_n(:,8);v_n(:,9);v_n(:,10)];
v_new=r+beta*Q*v;
%%%%%%%%%question 6
plot(a,pol_fn(1,:),a,pol_fn(2,:),a,pol_fn(3,:),a,pol_fn(4,:),a,pol_fn(5,:)),legend('z1','z2','z3','z4','z5');

%gini
p=[Mu(1,:);Mu(2,:);Mu(3,:);Mu(4,:);Mu(5,:)];
wealth=[a;a;a;a;a];
wg=gini(p,wealth,true);
title(['wealth gini index=',num2str(wg)]);

