function [values,dvalues] = BOOTerr(Data,B,Nrep,Name,Quantity,varargin)

% analysis of data by Bootstrap
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
% B        -- number of Bootsrap samples [100] (usual choices are ~25-200)
% Nrep     -- NOT USED (only for compatibility with UWerr, but I should add it)
% Name     -- NOT USED (only for compatibility with UWerr)
%------------------------------------------------
% output::
%  values  -- row vector with results (func applied to full sample)
% dvalues  -- their errors from fluctuation over jackknife-subsamples

[N,NQ]=size(Data);

if isempty(Quantity)  Quantity=1;    end
if isempty(B),        B=100;           end

% define a function that ...
if isnumeric(Quantity)
  func=@(d)d(:,Quantity); % ... works also when Quantity is a number ...
else
  func=@(d)Quantity(d,varargin{:}); % ... implicitely takes all the varargin.
end

%-------------------------------------------------------------------------

values=func(mean(Data,1));

for i =1:B
  BInd=ceil(N * rand(N,1)); % create a Boostrap index
  Fi(i,:) = func(mean(Data(BInd,:),1));
end

dvalues=sqrt(diag(cov(Fi)))';


