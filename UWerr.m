function [value,dvalue,ddvalue,tauint,dtauint,Qval] = ...
    UWerr(Data,Stau,Nrep,Name,Quantity,varargin)
%
% call:
% [value,dvalue,ddvalue,tauint,dtauint,Qval] = ...
% UWerr(Data,Stau,Nrep,Name,Quantity,P1,P2,....)
%
% autocorrelation-analysis of MC time-series
% following the Gamma-method in
% ``Monte Carlo errors with less errors''
% by Ulli Wolff, hep-lat/0306017
% 
%----------------------------------------------------------------------------
%  Ulli Wolff,   June 2004, Version 5
%  new: changes of hep-lat/0306017 v3 implemented concerning the error of tauint
%----------------------------------------------------------------------------
% input arguments beyond Data maybe omitted or set to []
% in the input list; they then get default value indicated in [D=..]
%
%
% Data     -- (N x Nalpha) matrix of measured (equilibrium!) data 
%              N = total number of measurements
%              Nalpha = number of (primary) observables
%              if your data are in a different format, consider 
%              the MATLAB commands `reshape' and `permute'
%
% Stau     -- guess for the ratio S of tau/tauint [D=1.5]
%             if set = 0, absence of autocorrelations is *assumed*
%
% Nrep     -- vector [N1 N2 N3 ...] specifying a breakup of the N rows
%             of Data into replica of length N1,N2,.. (N=sum(Nrep)!)  [D=[N]]
%             The replica distribution is histogrammed and a Q-value 
%             (=probability to find this much or more scatter) is given 
%             for R >= 2 replica; if you have one history you may set Nrep to 
%             artificially split it into >=2 replica to see their distribution                
%
% Name     -- if string: name of observable for titles of generated plots
%             if not string: all plots are supressed [D='NoName']
%
% Quantity -- either:
%             scalar function handle (@functionname) for the derived 
%             observable F; it has to operate on a row-vector of length 
%             Nalpha as first argument; optional parameters P1,P2,... are 
%             passed on to this function as 2nd, 3rd, .... argument
%          -- or:
%             integer between 1 and Nalpha to analyze primary observable [D=1]
%
%----------------------------------------------------------------------------
% output::
%  value   -- estimate for F
%  dvalue  -- its error (incl. autocorrelation effects)
%  ddvalue -- statistical error of the error
%  tauint  -- integrated autocorrelation time
%  dtauint -- statistical error of tauint
%  Qval    -- Q-value of the replica distribution if R >= 2
%             (goodness of fit to a constant)
%----------------------------------------------------------------------------
%

mlock  % this locks the function in memory and preserves persistent variables 
       % between calls; munlock('UWerr') required before a clear command 
persistent f1 f2 f3  % figure handles for plot windows
       
% analyze input arguments:

[N,Nalpha]=size(Data); % number of measurements and observables

if nargin < 5 | isempty(Quantity)  Quantity=1;    end
if nargin < 4 | isempty(Name),     Name='NoName'; end
if nargin < 3 | isempty(Nrep),     Nrep=[N];      end
if nargin < 2 | isempty(Stau),     Stau=1.5;      end

if any(Nrep ~= round(Nrep)) | any(Nrep < 1) | sum(Nrep) ~= N
  error('inconsistent N,Nrep')
end

if isnumeric(Quantity)
  if round(Quantity) ~= Quantity | Quantity < 1 | Quantity > Nalpha
    error('illegal numeric value of Quantity')
  end
  primary = 1; % logical variable
  primaryindex=Quantity;
else
  primary = 0;
end

Nrep=Nrep(:);   % make column
R=length(Nrep); % number of Replica

%----------------------------------------------------------------------------

% means of primaries:

abb=mean(Data,1);    %    total mean, primary obs.
abr=zeros(R,Nalpha); % replicum-mean, primary obs.
i0=1;
for r=1:R
  i1=i0-1+Nrep(r);
  abr(r,:)=mean(Data(i0:i1,:),1);
  i0=i0+Nrep(r);
end

% means of derived (or primary, depending on Quantity):

if primary
  Fbb=abb(primaryindex);
  Fbr=abr(:,primaryindex);
else
  Fbb=feval(Quantity,abb,varargin{:}); % total mean, derived obs.
  Fbr=zeros(R,1);                   % replica-means, derived obs.
  for i=1:R
    Fbr(i)=feval(Quantity,abr(i,:),varargin{:});
  end
end
Fb=sum(Fbr.*Nrep)/N;  % weighed mean of replica means

%----------------------------------------------------------------------------

% form the gradient of f and project fluctuations:

if primary
  delpro = Data(:,primaryindex)-abb(primaryindex);
else
  fgrad=zeros(Nalpha,1);
  h=std(Data,1)/sqrt(N); % spread for num. derivative
  ainc=abb;
  for alpha=1:Nalpha
    if h(alpha) == 0     % Data for this observable do not fluctuate
      fgrad(alpha)=0;
    else
      ainc(alpha) =abb(alpha)+h(alpha);
      fgrad(alpha)=feval(Quantity,ainc,varargin{:});
      ainc(alpha) =abb(alpha)-h(alpha);
      fgrad(alpha)=fgrad(alpha)-feval(Quantity,ainc,varargin{:});
      ainc(alpha) =abb(alpha);
      fgrad(alpha)=fgrad(alpha)/(2*h(alpha));
    end
  end
  
  % projected deviations: 
  delpro = Data*fgrad-abb*fgrad;
end

%----------------------------------------------------------------------------

% Computation of Gamma, automatic windowing
% note: W=1,2,3,.. in Matlab-vectorindex <-> W=0,1,2,.. in the paper!

% values for W=0:
GammaFbb(1) = mean(delpro.^2);
% sick case:
if GammaFbb(1) == 0
  beep
  fprintf('\n WARNING: no fluctuations \n')
  value   = Fbb;
  dvalue  = 0;
  ddvalue = 0;
  tauint  = 0.5;
  dtauint = 0;
  Qval=[];
  return
end

% compute Gamma, find the optimal W

if Stau == 0                    % no autocorrelations assumed
  Wopt=0;
  Wmax=0; flag=0;
else
  Wmax = round(min(Nrep)/2);    % do not take W larger than this
  flag=1; Gint=0;
end
W=1; 
while W <= Wmax,
  GammaFbb(W+1) = 0;
  % sum over replica:
  i0=1;
  for r=1:R
    i1=i0-1+Nrep(r);
    GammaFbb(W+1) = GammaFbb(W+1) + sum(delpro(i0:i1-W).*delpro(i0+W:i1));
    i0=i0+Nrep(r);
  end
  GammaFbb(W+1) = GammaFbb(W+1)/(N-R*W);
  if flag
    Gint=Gint+GammaFbb(W+1)/GammaFbb(1);
    if Gint <= 0, tauW=eps; else tauW = Stau/(log((Gint+1)/Gint)); end
    gW = exp(-W/tauW)-tauW/sqrt(W*N);
    if gW < 0                % this W is taken as optimal
      Wopt=W; 
      Wmax=min(Wmax,2*Wopt); flag = 0; % Gamma up to Wmax
    end
  end
  W=W+1;
end   % while-loop 

if flag,
  beep
  fprintf('\n WARNING: windowing condition failed up to W = %i \n',Wmax)
  Wopt=Wmax;
end

CFbbopt=GammaFbb(1) + 2*sum(GammaFbb(2:Wopt+1));   % first estimate
GammaFbb=GammaFbb+CFbbopt/N;                       % bias in Gamma corrected
CFbbopt=GammaFbb(1) + 2*sum(GammaFbb(2:Wopt+1));   % refined estimate
sigmaF=sqrt(CFbbopt/N);                            % error of F
tauintFbb=cumsum(GammaFbb)/GammaFbb(1)-0.5;

%----------------------------------------------------------------------------

% bias cancellation for the mean value

if R >= 2
  bF = (Fb-Fbb)/(R-1);
  Fbb=Fbb - bF;
  if abs(bF) > sigmaF/4
    beep
    fprintf('\n WARNING: a %.1f sigma bias of the mean has been cancelled \n',bF/sigmaF)
  end
  Fbr = Fbr - bF*N./Nrep;
  Fb  = Fb  - bF*R;
end

%----------------------------------------------------------------------------

% answers to be returned:

value   = Fbb;
dvalue  = sigmaF;
ddvalue = dvalue*sqrt((Wopt+0.5)/N);
tauint  = tauintFbb(Wopt+1);
% only change version 4 -> version 5 in next if-block:
if tauint > 0.5
  dtauint = tauint*2*sqrt((Wopt+0.5-1.5*tauint+0.125/tauint)/N);
else
  dtauint = tauint*2*sqrt((Wopt)/N);
end
% Q-value for replica distribution if R >= 2:

if R >= 2
  chisq=sum((Fbr-Fb).^2.*Nrep)/CFbbopt;
  Qval=1-gammainc(chisq/2,(R-1)/2);
else
  Qval=[];
end

%----------------------------------------------------------------------------

% Plotting:

% make plots, if demanded:

if ~ isstr(Name)   % no plots
  if ~ isempty(f1), close(f1), f1=[]; end
  if ~ isempty(f2), close(f2), f2=[]; end
  if ~ isempty(f3), close(f3), f3=[]; end
  return 
end

scrsz=get(0,'screensize'); % to place plots on the screen

% autocorrelation:

if Stau ~= 0
  
  if isempty(f1), 
    f1=figure; 
    set(f1,'Position',[0.01*scrsz(3) 0.45*scrsz(4) 0.45*scrsz(3) 0.45*scrsz(4)]);
    set(f1,'DefaultAxesFontSize',16);
  else 
    figure(f1)
  end; 
  clf;
  
  % plot of GammaF
  
  subplot(2,1,1)
  plot(0:Wmax,GammaFbb/GammaFbb(1),'x');                     % Gamma/Gamma(0)
  v=axis;v=[0 Wmax min(v(3),-0.1) 1];axis(v);hold on;
  plot([0 Wmax],[0 0],'--');                                 % zero line
  plot([Wopt Wopt],v(3:4),'r-');                             % Bar at Wopt
  title(['normalized autocorrelation of ' Name])
  ylabel('\Gamma')
  
  % plot of tauintF
  
  subplot(2,1,2)
  % tauint vs. W:
  errorbar(0:Wmax,tauintFbb,tauintFbb.*sqrt(([0:Wmax]+0.5-tauintFbb)/N)*2); 
  v=axis;v(1:2)=[0 Wmax];axis(v);hold on;
  plot([Wopt Wopt],v(3:4),'r-');               % Bar at Wopt
  plot([0 Wmax],[tauint tauint],'r--')         % hor. line at estimate
  title(['\tau_{int} with statistical errors of ' Name]);
  ylabel('\tau_{int}')
  xlabel('W')
  drawnow;
else
  close(f1); f1=[];
end

% histogram of data, if primary:

if primary
  if isempty(f2)
    f2=figure; 
    set(f2,'Position',[0.51*scrsz(3) 0.52*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4)]);
    set(f2,'DefaultAxesFontSize',16);
  else 
    figure(f2) 
  end
  clf
  hist(Data(:,primaryindex),20);
  title(['Distribution of all data for ' Name]);
  drawnow;
else
  close(f2); f2=[];
end

% histogram of replica values if R >= 2:

if R >= 2
  p=(Fbr-Fb)./(sigmaF*sqrt(N./Nrep-1));  
  if isempty(f3), f3=figure; else, figure(f3); end; clf;
  set(f3,'Position',[0.51*scrsz(3) 0.03*scrsz(4) 0.45*scrsz(3) 0.35*scrsz(4)]);
  set(f3,'DefaultAxesFontSize',16);
  hist(p);
  title(['replica distribution (mean=0,var=1) for ' Name ' >> Q = ' num2str(Qval,'%.2g') ' <<'])
  drawnow;
else
  close(f3); f3=[];
end

%----------------------------------------------------------------------------
%----------------------------------------------------------------------------
