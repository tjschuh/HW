function varargout=dummyplot(panelvars)
% [ah,ph,t,y]=DUMMYPLOT(panelvars)
%
% A simple function that plot four panels of stuff
%
% INPUT:
%
% panelvars
%
% OUTPUT:
%
% ah               Graphics handles to the four panels
% ph               Graphics handles to the four plot objects
% t                Some independent variable
% y                Some dependent variable
% 
% Originally written by tschuh@princeton.edu, 01/21/2021
% Last modified by tschuh@princeton.edu, 01/21/2021

% Defaults
defval('panelvars',[1 2 3 4])

% Calculate something
t=linspace(0,100,100);
for index=1:length(panelvars)
  y(:,index)=sin(2*pi*t(:)/panelvars(index));
end  

% Prepare figure
clf
cols={'k','b','g','r'};
lins={'-','--',':','-.'};
[ah,ha,H]=krijetem(subnum(2,2));

% Plot something
for index=1:length(panelvars)
  axes(ah(index))
  ph(index)=plot(t,y(:,index));
  xlabel('time')
end

% Cosmetics
for index=1:length(panelvars)
  axes(ah(index))
  ph(index).Color=cols{index};
  grid on
  longticks(ah(index))
  th(index)=title(sprintf('This is plot %i with parameter %3.1f',...
			  index,panelvars(index)));
  movev(th(index),0.05*range(ylim))
  nolabels(ah(1:2),1)
  nolabels(ha(3:4),2)
end

% Optional outputs
varns={ah};
varargout=varns(1:nargout);
