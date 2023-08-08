function handle = imab(data, clims, bSep, noCols)
% function handle = imab(data, clims, bSep, noCols)
%
% Just show the 2D, 3D or 4D data as quickly as possible
% 
% Concatenates the slice dimension from left to right, and the time
% dimension from top to bottom.
%
% When Mosaiced, lowest slice is in bottom-left corner.
%
% bSep is boolean option to separate each image from its neighbours rather
% than making one big image (when Mosaicing)
%
% noCols is no. of columns to use in Mosaic (only used for single volume)
%
% January 2010
% daniel.gallichan@uniklinik-freiburg.de
%
% subfunction 'makemosaic' borrowed from Jaemin Shin's implementation.


if nargin < 4; noCols =[]; end
if nargin < 3; bSep = 0; end
if nargin < 2
    clims = [];
end

if size(data,3)==1
    data = squeeze(data);
end

if size(data,4) > 20
    disp('Dataset more than 20 time points, could be difficult to see what''s going on...!')
end

data = permute(data,[2 1 3 4]);
if ~isreal(data)
    data = abs(data);
    disp('Data were complex - absolute values are shown')
end



if bSep && size(data,3)>1
    figure
    dims = size(data);
    
    if isempty(clims)
        clims = [min(data(:)) max(data(:))];
    end
        
    if size(data,4)==1
        if isempty(noCols)
            noCols = ceil(dims(3).^0.5);%min(5,ceil(dims(3)/2));
        end
        noRows = ceil(dims(3)/noCols);
        iPlots = reshape(noCols*noRows:-1:1,[noCols noRows]).';
        iPlots = iPlots(:,end:-1:1);
        
  
        
        subplot1(noRows,noCols)
        for iR = 1:noRows
            for iC = 1:noCols
                if (noRows==1), 
                    subplot1(iC)
                else
                    subplot1([iR iC])
                end
                if iPlots(iR,iC) > dims(3)
                    set(gca,'visible','off')
                else
                    hndl=imagesc(data(:,:,iPlots(iR,iC)),clims);axis equal tight;set(gca,'Xtick',[],'YTick',[]);
                end
            end
        end
        

    else

        noRows = size(data,4); noCols = size(data,3);
        subplot1(noRows,noCols)
        for iR = 1:noRows
            for iC = 1:noCols
                if (noRows==1),
                    subplot1(iC)
                else
                    subplot1([iR iC])
                end
                hndl=imagesc(data(:,:,iC,iR),clims);axis equal tight;set(gca,'Xtick',[],'YTick',[]);
            end
        end

    end

else
    %%%% Normal behaviour:-
    if size(data,4)==1
        if size(data,3)>1
            if isempty(noCols)
                noCols = ceil(size(data,3).^0.5);%min(5,ceil(size(data,3)/2));
            end
            outIm = makemosaic(data,noCols);
        else
            outIm = data;
        end
        if isempty(clims)
            hndl=imagesc(outIm);axis equal tight;set(gca,'Xtick',[],'YTick',[]);
        else
            hndl=imagesc(outIm,clims);axis equal tight;set(gca,'Xtick',[],'YTick',[]);
        end
    else
        if isempty(clims)
            hndl=imagesc(imcat(imcat(data,2,3),1,4,0));axis equal tight;set(gca,'Xtick',[],'YTick',[]);
        else
            hndl=imagesc(imcat(imcat(data,2,3),1,4,0),clims);axis equal tight;set(gca,'Xtick',[],'YTick',[]);
        end
    end

end

set(gca,'YDir','normal');

if nargout
    handle=hndl;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function imall=makemosaic(im,MaxN);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function imall=makemosaic(im,MaxN);
% Make a mosaic image for display the image with "show.m"
% i.e., the 3D image [im] transforms to a mosaic 2D image [imall]
% If [im] is 4D, [im(:,:,:,1)] will be used
% NOTE : First, singleton dimensions will be removed; 64x64x1x20 -> 3D
% Input :
%   [im] : 3D or 4D image
%   [MaxN](option): The number of colons, Default is 5
% Output :
%   [imall]: mosaic 2D image
% Usages,
% imall=makemosaic(im,MaxN);
% imall=makemosaic(im);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Examples,
% let [im] is a matrix of 64x64x20
% imall=makemosaic(im,10);
% [imall] is a 2x10 image of 64x64, size(imall)= 128 x 640
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright (c) Jaemin Shin
% jaemins@gatech.edu
% 01/30/06
% updated 
% 02/28/06 : bug fixed

if exist('MaxN','var') == 0
    MaxN = 5;
end
im = squeeze(im);
dim = size(im);
if length(dim) < 2;
    error('Input is 1D or 2D signal')
elseif length(dim) ==4
    im = squeeze(im(:,:,:,1));
    disp('4D : TimePoint 1 was used')
elseif length(dim) > 4
    error('5D or Higher dimension does not support')
end
Nrow = ceil(dim(3)/MaxN);
Rcol = mod(MaxN - mod(dim(3),MaxN),MaxN);

if dim(3) <= MaxN
    imall = reshape(im,[dim(1) dim(2)*dim(3)]);
    imall = [imall,zeros(dim(1),dim(2)*Rcol)];
else
    imall = reshape(im(:,:,1:MaxN),[dim(1) dim(2)*MaxN]);
    for ii=2:Nrow-1 % bug fixed
        temp = reshape(im(:,:,(ii-1)*MaxN+1:ii*MaxN),[dim(1) dim(2)*MaxN]);
        imall = cat(1,imall,temp);
    end
    temp = reshape(im(:,:,(Nrow-1)*MaxN+1:end),[dim(1) dim(2)*(MaxN-Rcol)]);
    temp = [temp,zeros(dim(1),dim(2)*Rcol)];
    imall = cat(1,imall,temp);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function out = imcat(data, outdim, indim, updown)
% function out = imcat(data, outdim, indim, updown)
% concatenate data along indim dimension and add it to dimension outdim

if nargin < 4  updown = 1;    end
if nargin < 3  indim = 3;     end

if updown
    switch indim
        case 3
            temp = data(:,:,1,:);
            if size(data,3)> 1
                for i = 2:size(data,3)
                    temp = cat(outdim,temp, data(:,:,i,:));
                end
            end

        case 4
            temp = data(:,:,:,1);
            if size(data,4)> 1
                for i = 2:size(data,4)
                    temp = cat(outdim,temp,data(:,:,:,i));
                end
            end

        case 5
            temp = data(:,:,:,:,1);
            if size(data,5)> 1
                for i = 2:size(data,5)
                    temp = cat(outdim,temp,data(:,:,:,:,i));
                end
            end

    end

else
    switch indim
        case 3
            temp = data(:,:,end,:);
            if size(data,3)> 1
                for i = size(data,3)-1:-1:1
                    temp = cat(outdim,temp, data(:,:,i,:));
                end
            end

        case 4
            temp = data(:,:,:,end);
            if size(data,4)> 1
                for i = size(data,4)-1:-1:1
                    temp = cat(outdim,temp,data(:,:,:,i));
                end
            end

        case 5
            temp = data(:,:,:,:,end);
            if size(data,5)> 1
                for i = size(data,5)-1:-1:1
                    temp = cat(outdim,temp,data(:,:,:,:,i));
                end
            end

    end
end

axis off
set(gca,'box','off')
out = temp;