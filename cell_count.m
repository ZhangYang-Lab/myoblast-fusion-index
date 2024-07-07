clear 
clc

SmallestMyotubePixelCount = 500;
SmallestNucleusPixelCount = 300;
MyotubeChannel = 1;
FillSize = 10;
NucFillSize = 5;
TubeThresh = 1.0;
MaxNucSizeDivisor = 100;

numerics = ["SmallestMyotubePixelCount","SmallestNucleusPixelCount",...
    "MyotubeChannel","NucChannel","FillSize","NucFillSize","MinCircleRad",...
    "MaxCircleRad","TubeThresh","MinNuclei","MaxNucSizeDivisor"];

[FileName,PathName] = uigetfile('*.tif','Please select tif');
FilePath = fullfile(PathName, FileName);

%if there's a wildcard in the path given, find all files matching it and
%run tool against all sequentially.
fileList = dir(FilePath);
filecount = size(fileList);
disp(fileList);

File = fullfile(fileList(1).folder, fileList(1).name);
disp(File);

img = imread(File); % Read image
imshow(img);

tubes = getmyotubes(img,MyotubeChannel,SmallestMyotubePixelCount,FillSize,TubeThresh);
nuclei = getnuclei(img,MyotubeChannel,SmallestNucleusPixelCount,NucFillSize,MaxNucSizeDivisor);

mask = imbinarize(tubes + nuclei);

filled = imfill(mask,'holes');
opened = imopen(filled, ones(FillSize,FillSize));
%ignore anything smaller than 500 pixels.
openmask = bwareaopen(opened, SmallestMyotubePixelCount);
closemask = imclose(openmask,strel('disk',30));
filled = imfill(closemask,'holes');

[x,y]=size(closemask);
pixels=x*y;
mask = bwareafilt(filled,[SmallestMyotubePixelCount pixels]);
Finalmask(:,:,1) = mask;
Finalmask(:,:,2) = mask;
Finalmask(:,:,3) = mask;
Finalmask = uint8(Finalmask);

figure(2)
imshow(mask)

[FileName,PathName] = uigetfile('*.tif','Please select tif');
FilePath = fullfile(PathName, FileName);
nuclei_img = imread(FilePath); % Read image
figure(3)
imshow(nuclei_img);

nuclei_in_tube = nuclei_img .* Finalmask;
figure(4)
imshow(nuclei_in_tube)

figure(5)
imshow(nuclei_in_tube+img)

nuclei_counts_full = getnucleicounts(nuclei_img, SmallestMyotubePixelCount, TubeThresh, FillSize);
nuclei_counts_in_tube = getnucleicounts(nuclei_in_tube, SmallestMyotubePixelCount, TubeThresh, FillSize);
ratio = nuclei_counts_in_tube / nuclei_counts_full;

disp_info = ['nuclei_counts_full is :',num2str(nuclei_counts_full), ...
    ' nuclei_counts_in_tube is : ', num2str(nuclei_counts_in_tube), ...
    ' ratio is : ', num2str(ratio),];

disp(disp_info)
%%
function y = getmyotubes(image,MyotubeChannel,SmallestMyotubePixelCount,FillSize,thresh)
    singlecolour = image(:,:,MyotubeChannel);
    
    I_cropped = singlecolour;
    I_eq = adapthisteq(I_cropped);
    
    % mark myotubes
    marked = imbinarize(I_eq, graythresh(I_eq)*thresh);
    filled = imfill(marked,'holes');
    opened = imopen(filled, ones(FillSize,FillSize));
    %ignore anything smaller than 500 pixels.
    openmask = bwareaopen(opened, SmallestMyotubePixelCount);
    
    [x,y]=size(I_eq);
    pixels=x*y;
    
    
    Finalmask = bwareafilt(openmask,[SmallestMyotubePixelCount pixels]);
    
    %Split into parts
    [indeximgarray,~] = bwlabel(Finalmask);
    
    y = indeximgarray;
end

% 
function y = getnuclei(image,NucChannel,SmallestNucleusPixelCount,NucFillSize,MaxNucSizeDivisor)
    %get the correct channel
    singlecolour = image(:,:,NucChannel);
    I_eq = adapthisteq(singlecolour);
    I_eq = wiener2(I_eq,[NucFillSize NucFillSize]);
    
    % use multiple thresholds, extract an image of the scene for each threshold
    
    thresh = multithresh(I_eq,2);
    seg_I = imquantize(I_eq,thresh);
    map = [1, 0, 0
        0, 1, 0
        0, 0, 1];
    RGB = label2rgb(seg_I,map);
    
    
    singlecolour1 = ~logical(RGB(:,:,1));
    singlecolour2 = ~logical(RGB(:,:,2));
    singlecolour3 = logical(RGB(:,:,3));
    
    %extract only elements which are small enough to be nuclei
    
    [x,y]=size(I_eq);
    pixels=x*y;
    %throw out anything larger than 1% of the total image size to discard
    %washed out regions
    upper = pixels/MaxNucSizeDivisor;
    %throw out anything less than this many pixels
    lower = SmallestNucleusPixelCount;
    
    BW2 = bwareafilt(singlecolour1,[lower upper]);
    BW3 = bwareafilt(singlecolour2,[lower upper]);
    BW4 = bwareafilt(singlecolour3,[lower upper]);
    
    [indeximgarray,~] = bwlabel(BW2 + BW3 + BW4);
   
    
    y = indeximgarray;
end

%
function y = splitconnected(shapes)
    %credit https://blogs.mathworks.com/steve/2013/11/19/watershed-transform-question-from-tech-support/
    D = -bwdist(~shapes);
    mask = imextendedmin(D,2);
    D2 = imimposemin(D,mask);
    Ld2 = watershed(D2);
    test2 = shapes;
    test2(Ld2 == 0) = 0;
    y = test2;
end

%
function num_cells = getnucleicounts(img, SmallestMyotubePixelCount, thresh, FillSize)

% 转换为灰度图像
gray_img = rgb2gray(img);

I_eq = adapthisteq(gray_img);
marked = imbinarize(I_eq, graythresh(I_eq)*thresh);
filled = imfill(marked,'holes');
opened = imopen(filled, ones(FillSize,FillSize));
%ignore anything smaller than 500 pixels.
openmask = bwareaopen(opened, SmallestMyotubePixelCount);

[x,y]=size(I_eq);
pixels=x*y;


Finalmask = bwareafilt(openmask,[SmallestMyotubePixelCount pixels]);
bw4 = splitconnected(Finalmask);

[label_img, num_cells] = bwlabel(bw4);

% 获取每个细胞核的中心位置
stats = regionprops('table', label_img, 'Centroid');
centroids = stats.Centroid;

% 在原始图像中显示细胞核位置
figure;
imshow(img);
hold on;
plot(centroids(:,1), centroids(:,2), 'r*', 'MarkerSize', 10);
title(['分割后的细胞核数量：', num2str(num_cells)]);

end


