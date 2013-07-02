function [values,dvalues] = JACKerr(Data,B,Nrep,Name,Quantity,varargin)

% analysis of MC time-series by jackknife binning
%------------------------------------------------
% input::
% Data     -- matrix, where consecutive ROWS of length m correspond 
%              to consecutive measurements of m observables
% Quantity -- either:
%             function handle (@functionname) for the derived 
%             observable F; it has to operate on a row-vector of length 
%             NQ as first argument; optional parameters P1,P2,... are 
%             passed on to this function as 2nd, 3rd, .... argument
%             the output is of the form [1,nv]=size(Quantity(Data,P1,P2,...))
%          -- or:
%             integer between 1 and NQ to analyze primary observable [D=1]
%
% B        -- number of measurements combined to bins 
% Nrep     -- NOT USED (only for compatibility with UWerr, but I should add it)
% Name     -- NOT USED (only for compatibility with UWerr)
%------------------------------------------------
% output::
%  values  -- row vector with results (func applied to full sample)
% dvalues  -- their errors from fluctuation over jackknife-subsamples

[N,NQ]=size(Data);

if isempty(Quantity)  Quantity=1;    end
if isempty(B),        B=1;           end

% define a function that ...
if isnumeric(Quantity)
  func=@(d)d(:,Quantity); % ... works also when Quantity is a number ...
else
  func=@(d)Quantity(d,varargin{:}); % ... implicitely takes all the varargin.
end

%-------------------------------------------------------------------------

Nb=floor(N/B); % number of bins.
n0 = N-Nb*B; % discard the first n0 measures that do not fit the binning choice.
Data1=Data(n0+1:N,:);

values=func(mean(Data1,1));

if(Nb>1)
  for i =1:Nb
    bin=[1+(i-1)*B:i*B];
    cbin=setdiff([1:N-n0],bin);
    Fi(i,:) = func(mean(Data1(cbin,:),1));
  end
dvalues=sqrt(diag(cov(Fi))'/Nb)*(Nb-1);
end

