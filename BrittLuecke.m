function [dev dz chi2 dxG dzG W A] = BrittLuecke(x,t,z,sigma,F,out_flag,flag)
% error analysis following: 
%"The Estimation of Parameters in Nonlinear Implicit Models"
% By H.I. Britt and R.H. Luecke (1973)
% Matlab code by Luigi Scorzato (2007)
% x:     measured values for the Q "regressors"(e.g. m_pcac,om,Z's) and responses (m_pi,f_pi..))
%        in  N simulation points. [Qr rows, N cols].
% t:     variables of F with no errors.
% [M,N] = size out F
% sigma: N covariance matrices [Q rows,Q cols each].
% z:     parameters (B_0, f, LEC's, W's) [P rows, 1 cols].
% F:     handel to functions (laws) of x and z, which  must accept
%        col vector z of size [P,1], matrices x of size [Q,N]
%        and return ~0 of size [M,N].
% flag = 1: compute ordinary CHI2 (BrittLuecke Eq. (2)) [dev]
% flag = 2: compute Generalized square deviation        [dev]
% flag = 3: compute Generalized square deviation
%           and parameters  errors                      [gdev,dz].
% flag > 3: compute errors [dz]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sizes and initialization.
ep=6e-6;
P=size(z,1);
[T,N] =size(t);
[Q,N]=size(x);
v=feval(F,x,t,z);
[M,N]=size(v);
dev=0;
ldev=0;
dz=zeros(P,1);
dzG=zeros(M,P,N);
dxG=zeros(M,Q,N);
W=zeros(M,M,N);
A=zeros(P,P);
persistent constshift

% naive weighted square deviation
% assumes that the output of F corresponds to the first M rows-columns of sigma.
for j=1:N;
  isigma(:,:,j)=inv(sigma(:,:,j));
  s(j)= v(:,j)' * isigma(find(out_flag),find(out_flag),j) * v(:,j);
end
chi2(1)=sum(s);

% compute dxG
tt=reshape(t,T,1,N);
xx=reshape(x,Q,1,N);
h=xx.*ep;
h(find(xx==0))=ep;
for j=1:N;
  tl = tt(:,:,j)*ones(1,Q);
  xp = xx(:,:,j)*ones(1,Q) + diag(h(:,:,j));
  xm = xx(:,:,j)*ones(1,Q) - diag(h(:,:,j));
  dh = ones(M,1)*h(:,:,j)';
  dxG(1:M,1:Q,j)= -(feval(F,xp,tl,z)-feval(F,xm,tl,z))./(2.*dh);
end
%compute Gen square dev.
for j=1:N
  tempW=(dxG(:,:,j) * sigma(:,:,j) * dxG(:,:,j)');
  W(:,:,j)=inv(tempW);
  gs(j)=v(:,j)' * W(:,:,j) * v(:,j);
  ldw(j)=-log(det(W(:,:,j)));
end
ldev= sum(ldw);
if(isempty(constshift))
  constshift=-ldev;
end
chi2(2)=sum(gs);
chi2(3)=chi2(2)+ldev+constshift; 
NMP=max(N*M-P,1);
if(flag<=3)
  dev=chi2(flag)/NMP;
else
  dev=chi2(1)/NMP;
end
[dev NMP]

if(flag>3)
  % compute dzG
  k=z.*ep;
  k(find(z==0))=ep;
  zp = z*ones(1,P) + diag(k);
  zm = z*ones(1,P) - diag(k);
  for j=1:P;
    dk = ones(M,N)*k(j); 
    dzGt(:,:,j)= -(feval(F,x,t,zp(:,j))-feval(F,x,t,zm(:,j)))./(2.*dk);
  end
  dzG=permute(dzGt,[1,3,2]);
  % compute dz
  for j=1:N
    A(:,:)=A(:,:)+(dzG(:,:,j)' * W(:,:,j) * dzG(:,:,j));
  end
  dz=sqrt(diag(inv(A)));
end
