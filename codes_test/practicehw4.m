clear;
close all;
%%%%set up parameters
b=0.994;
s=1.5;
pi=[0.97, 1-0.97;1-0.5,0.5];
pi_inv=pi^100;
b=0.5;
y=[1, b];

%%%%set up grids
a_min=-2;
a_max=5;
a_num=1000;
a=linspace(a_min,a_max,a_num);


%%%consumption and return function

q_min=0;
q_max=1;

q_guess=(q_min+q_max)/2;
cons=bsxfun(@minus,a',q_guess*a);
cons=bsxfun(@plus,cons,permute(y,[1 3 2]));
ret=(cons.^(1-s))/(1-s);
ret(cons<0)=-Inf;

%%%value function iteration
dis=1;
tol=1e-6;
while dis>tol
    v_guess=zeros(2,a_num);
    vf=ret+beta*repmat(permute(pi*v_guess,[3 2 1]),[a_num 1]);
    [vfn, policy_ind]=max(vf,(),2];
    vfn=shiftdim(vfn,2];
    dis=max(abs(vfn(:)-v_guess(:)));
    v_guess=vfn;
end


policy_ind=shiftdim(policy_ind,2);

%%%%%%%distribution iteration
u=ones(2,a_num)/(2*a_num);
[emp_ind,a_ind]=find(u>0);
dif=1;
while dif>1e-6
    u_new=
