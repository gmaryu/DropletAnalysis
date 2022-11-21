%function overlayTrackingResult(sumTrackTable,PosList,genPos)
%{
    
    Generate overlay image of label image and original image
    INPUTS
        path_label: path to label folder
        path_img: path to image folder
        path_save: path to ovelay image
        tracks: target position's sumTrackTable
        frameNum: the number of timepoints
    Segmentation results and tracking results mat files must be loaded.
%}
tracks=sumTrackTable{1};
frameNum=500;

% label image
label_folder='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\';
path_label=[label_folder,'20220128_CycB-mNG/SegmentationResults/Pos0_segmentation/'];
%nucFormat='img_000000%03d_nuc_000.jpg';
maskFormat='img_000000%03d_segment_000.jpg';
disp(path_label);

% original image
image_folder='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\data\';
path_img=[image_folder,'20220128_CycB-mNG/Pos0/'];
imgFormat='img_000000%03d_5-GFP_000.tif';
disp(path_img);

% move to save path
label_folder='\\biop-qiongy-nas.biop.lsa.umich.edu\qiongy-data\users\Gembu\results\';
path_save=[label_folder,'20220128_CycB-mNG/TrackingResults/Pos0/'];
if ~exist(path_save, 'dir')
       mkdir(path_save)
end
cd(path_save);
cellids=extractfield(tracks, 'id');
for i=1:frameNum
    %disp(i);
    %% make overlay image
    ifname= sprintf(imgFormat,i-1);
    mfname= sprintf(maskFormat,i-1);
    full_imgPath=[path_img,ifname];
    full_nucPath=[path_label,mfname];
    img= imread(full_imgPath);
    img= double(img);
    img=img / max(img(:)) * 255;
    %img= img./max(img(:))*255;
    nucimg= imread(full_nucPath);
    se = strel('disk',2);
    nucimg=imopen(nucimg, se);
    olimg=labeloverlay(uint8(img), logical(nucimg));
    
    %% add tracking ID to overlay image
    text_position=[];
    text_str={};
    box_color={};
    for j=1:length(cellids)
        cellid=cellids(j); %actual droplet id
        track=tracks(j);
        if sum(track.Feat.t==i)==1
            %disp(j);
            tmp_index=find(track.Feat.t==i);
            tmp_xcoord=round(track.Feat.xcoord(tmp_index))-10;
            tmp_ycoord=round(track.Feat.ycoord(tmp_index))-10;
            text_position=vertcat(text_position, [tmp_xcoord,tmp_ycoord]);
            text_str=horzcat(text_str,num2str(j));
            box_color = horzcat(box_color,'yellow');
        end
    end
    olimg = insertText(olimg,text_position,text_str,'FontSize',18,'BoxColor',...
    box_color,'BoxOpacity',0.1,'TextColor','white');
    %% save image
    nameFormat1='img_000000%03d_Tracking_Overlay_000.jpg';
    fname1= sprintf(nameFormat1,i-1);
    imwrite(olimg, fname1);
end
