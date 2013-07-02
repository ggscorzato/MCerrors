function errorbarxy(x,y,lx,ly,ux,uy,linecol,errorcol,linesty,errorsty,linewid,errorwid)

%errorbarxy(x,y,lx,ly,ux,uy,linecol,errorcol,linesty,errorsty,linewid,errorwid)
%James Rooney,  17 October 2003
% modified by Luigi Scorzato, 16 August 2005

%This function allows the user to plot the graph of x against y,
%along with both x and y errorbars.
%For the x and y errors it is possible to input both lower (lx and
%ly)  and upper  (ux and uy) values for the
%errors at a particular point.  If the upper values are not
%specified then the program assumes the errors 
%are symmetrical and use the lower values.  it is also possible to
%specify the plot line colour, marker, and 
%linestyle using the standard 'plot' command notation in the input
%variable 'linecol'.  Also the line colour for 
%the errobars can be specified in the variable 'errorcol'.  It is
%important to note that if these colour options 
%are to be used and any of the error limit vectors are empty then
%they should not be excluded, but presented 
%in a [] form signifying an empty vector.
%

if exist('linecol','var')==0 | isempty(linecol)
  linecol='b';
end

if exist('linesty','var')==0 | isempty(linesty)
  linesty='.';
end

if exist('linewid','var')==0 | isempty(linewid)
  linewid=0.5;
end

if exist('errorcol','var')==0 | isempty(errorcol)
  errorcol='r';
end

if exist('errorsty','var')==0 | isempty(errorsty)
  errorsty='-';
end

if exist('errorwid','var')==0 | isempty(errorwid)
  errorwid=0.5;
end

plot(x,y,'color',linecol,'LineStyle',linesty,'LineWidth',linewid,'MarkerSize',16)
hold on

xw=(max(x)-min(x))/100;
yw=(max(y)-min(y))/100;


lye=exist('ly','var');
lxe=exist('lx','var');
uye=exist('uy','var');
uxe=exist('ux','var');

if lye+lxe+uye+uxe==0 | isempty(lx) & isempty(ux) & isempty(ly) & isempty(uy)
  return
end

if uye==0 | isempty(uy)
  uy=ly;
end

if uxe==0 | isempty(ux)
  ux=lx;
end

for t=1:length(x)
  
  if ~isempty(ux)
    %x errorbars
    line([x(t)-lx(t) x(t)+ux(t)],[y(t) y(t)],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
    line([x(t)-lx(t) x(t)-lx(t)],[y(t)-yw y(t)+yw],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
    line([x(t)+ux(t) x(t)+ux(t)],[y(t)-yw y(t)+yw],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
  end
  if ~isempty(uy)
    %y errorbars
    line([x(t) x(t)],[y(t)-ly(t) y(t)+uy(t)],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
    line([x(t)-xw x(t)+xw],[y(t)-ly(t) y(t)-ly(t)],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
    line([x(t)-xw x(t)+xw],[y(t)+uy(t) y(t)+uy(t)],'color',errorcol, ...
	 'LineStyle',errorsty,'LineWidth',errorwid)
  end    
end
hold off
    
