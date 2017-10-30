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
a_num=500;
a=linspace(a_min,a_max,a_num);
%%%%%%question 5%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nd=ns;
k_guess=2;

d=1;
while d>=0.001;
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
v_tol = 1;
    while v_tol >.0001;
        % CONSTRUCT RETURN + EXPECTED CONTINUATION VALUE
      
        vf=bsxfun(@plus,ret,permute(beta*zprob*v_guess,[3,2,1]));
        
        % CHOOSE HIGHEST VALUE (ASSOCIATED WITH a' CHOICE)
        [vfn,pol_indx]=max(vf,[],2);
        v_tol=[max(abs(vfn(:,:,1)' - v_guess(1,:))) ; max(abs(vfn(:,:,2)' - v_guess(2,:)));...
            max(abs(vfn(:,:,3)' - v_guess(3,:)));max(abs(vfn(:,:,4)' - v_guess(4,:)));...
            max(abs(vfn(:,:,5)' - v_guess(5,:)))];
        v_guess=[vfn(:,:,1)';vfn(:,:,2)';vfn(:,:,3)';vfn(:,:,4)';vfn(:,:,5)'];
    end;
  % KEEP DECSISION RULE
    pol_indx=permute(pol_indx, [3 1 2]);
    pol_fn=a(pol_indx); 
    % SET UP INITITAL DISTRIBUTION
    Mu=ones(m,a_num)/(m*a_num);
    
    % ITERATE OVER DISTRIBUTIONS
     m_tol=1;
     while m_tol>0.0001
        [z_ind, a_ind, mass] = find(Mu > 0); % find non-zero indices
    
        MuNew = zeros(size(Mu));
   
        
       for ii = 1:length(z_ind)
        apr_ind = pol_indx(z_ind(ii), a_ind(ii)); % which a prime does the policy fn prescribe?
        MuNew(:, apr_ind) = MuNew(:, apr_ind) +[zprob(z_ind(ii),1)*Mu(z_ind(ii),a_ind(ii));zprob(z_ind(ii),2)*Mu(z_ind(ii),a_ind(ii));zprob(z_ind(ii),3)*Mu(z_ind(ii),a_ind(ii));zprob(z_ind(ii),4)*Mu(z_ind(ii),a_ind(ii));zprob(z_ind(ii),5)*Mu(z_ind(ii),a_ind(ii))];
        % which mass of households goes to which exogenous state?
       end
       m_tol=max(max(abs(MuNew-Mu)));
       Mu=MuNew;
     end     

 aggsav=Mu(1,:)*a'+Mu(2,:)*a'+Mu(3,:)*a'+Mu(4,:)*a'+Mu(5,:)*a';
 d=abs(k_guess-aggsav);
 if k_guess<aggsav
     k_guess=k_guess+abs(k_guess-aggsav)/2;
 else k_guess=k_guess-abs(k_guess-aggsav)/2;
 end
 k_guess;
 d;
end

%%%%%%%%%question 6
plot(a,pol_fn(1,:),a,pol_fn(2,:),a,pol_fn(3,:),a,pol_fn(4,:),a,pol_fn(5,:)),legend('z1','z2','z3','z4','z5');

%gini
p=[Mu(1,:);Mu(2,:);Mu(3,:);Mu(4,:);Mu(5,:)];
wealth=[a;a;a;a;a];
wg=gini(p,wealth,true);
title(['wealth gini index=',num2str(wg)]);

