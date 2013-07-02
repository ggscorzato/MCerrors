function [dev chi2 dz W dzG A] = MBerr(y,x,t,z,sigma,F,flag)
% Error analysis following: 
%  "Generalized least-squares fit of multiequation models" S.L.Marshall and J.G.Blencoe.
% Author: Luigi Scorzato (2007)
% x:     measured values for the Qr "regressors"(e.g. m_pcac,om,Z's)
%        in  N simulation points. [Qr rows, N cols].
% t:     variables of F with no errors.
% y:     measured values for the M=Q-Qr responses (e.g. m_pi,f_pi,g_pi)
% sigma: N covariance matrices [Q=M+Qr rows,Q=M+Qr  cols each].
% z:     parameters (B_0, f, LEC's, W's) [P rows, 1 cols].
% F:     handel to functions (laws) of x and z, which  must accept
%        col vector z of size [P,1], matrices x of size [Qr,*] 
%        and return y1 of size [M,*].
% flag = 0: compute naive (weighted) square deviation  [dev]
% flag = 1: compute Generalized square deviation       [dev]
% flag = 2: compute Generalized square deviation
%           and parameters  errors                     [gdev,dz].
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% sizes and initialization.
ep=6e-6;
P=size(z,1);
[T,N] =size(t);
[Qr,N]=size(x);
y1=feval(F,x,t,z);
g=y-y1;
[M,N]=size(y1);
Q=Qr+M;
fullx=[y;x]; %(although not used, it clarifies the order.)
dev=0;
ldev=0;
dz=zeros(P,1);
W=zeros(M,M,N);
dzG=zeros(M,P,N);
A=zeros(P,P);
persistent constshift
%constshift=10*N*M; % This formula has no deep meaning: I simply need a constant as close as possible
                     % to log(det): NM make sense, 8 ... no idea.

% naive weighted square deviation
v=(y-y1);
for j=1:N;
  isigma(:,:,j)=inv(sigma(:,:,j));
  s(j)= v(:,j)' * isigma(1:M,1:M,j) * v(:,j);
end
NMP=N*M-P;
if(NMP>0)
  dev=sum(s) / NMP;
else
  dev=sum(s);
end
if(flag>0)
  % compute dxG
  tt=reshape(t,T,1,N);
  xx=reshape(x,Qr,1,N);
  h=xx.*ep;
  h(find(xx==0))=ep;
  for j=1:N;
    dxG(1:M,1:M,j)=eye(M,M);
    tl = tt(:,:,j)*ones(1,Qr);
    xp = xx(:,:,j)*ones(1,Qr) + diag(h(:,:,j));
    xm = xx(:,:,j)*ones(1,Qr) - diag(h(:,:,j));
    dh = ones(M,1)*h(:,:,j)';
    dxG(1:M,M+1:Q,j)= -(feval(F,xp,tl,z)-feval(F,xm,tl,z))./(2.*dh);
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
  chi2=sum(gs);
  dev=chi2+ldev+constshift;
  if(NMP>0)
    dev=dev/NMP;
    chi2=chi2/NMP;
  else
    dev=dev;
    chi2=chi2;
  end
end
%[dev chi2]
[dev]

if(flag>1)
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
%  dz=diag(inv(A)) *sqrt(chi2); WRONG in Marshall Blencoe.
  dz=sqrt(diag(inv(A)));
end
