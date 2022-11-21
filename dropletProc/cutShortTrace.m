function tracks1=cutShortTrace(tracks,frameNum,thres)
Len = [tracks.Len].';
tracks1=tracks(Len>(frameNum*thres));

for i = 1:size(tracks1, 1)
    %disp(i)
    tracks1(i).id=i;
end
end