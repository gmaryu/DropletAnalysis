function S1=skwDroplet(S1)
for jj=1:size(S1)
    S1(jj).skewness=skewness(double(S1(jj).PixelValues));
    S1(jj).stdInt=std(double(S1(jj).PixelValues));
end
end
