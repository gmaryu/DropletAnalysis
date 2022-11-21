function [lineOut, fillOut] = stdshade(amatrix,alpha,acolor,F)
%% stdshade
% make a shade of std value over the mean value
% usage: stdshading(amatrix,alpha,acolor,F,smth)

% - acolor defines the used color (default is red) 
% - F assignes the used x axis (default is steps of 1).
% - alpha defines transparency of the shading (default is no shading and black mean line)

if exist('acolor','var')==0 || isempty(acolor)
    acolor='r'; 
end

if exist('F','var')==0 || isempty(F)
    F=1:size(amatrix,2);
end


if ne(size(F,1),1)
    F=F';
end

amean = nanmean(amatrix,1); %get man over first dimension
astd = nanstd(amatrix,[],1); % to get std shading

if exist('alpha','var')==0 || isempty(alpha) 
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor,'linestyle','none');
    acolor='k';
else
    fillOut = fill([F fliplr(F)],[amean+astd fliplr(amean-astd)],acolor, 'FaceAlpha', alpha,'linestyle','none');
end

if ishold==0
    check=true; else check=false;
end

hold on;
lineOut = plot(F,amean, 'color', acolor,'linewidth',1.5); %% change color or linewidth to adjust mean line

if check
    hold off;
end

end