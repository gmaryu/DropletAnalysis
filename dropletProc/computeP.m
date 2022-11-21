function p=computeP(x,delta)
% major probability model function for feature
%parWeight=0.1;
% p1=exp(-x(1).*x(1)./delta(1).^2);
% p=parWeight+(1-parWeight)*prod(exp(-x(2:end).^x(2:end)./delta(2:end).^2));
% p=p*p1;
p=prod(exp(-x.*x./delta.^2));