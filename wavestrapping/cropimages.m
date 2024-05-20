%% Input
imName = {'1';...
    '2';'3';'4';'5';'6';...
    '7'; '8';'9';'10';'11';'12';...
    '13';'14';'15'};

%% Input
for j = 1:15


inImage = imread(strcat('stimuli_color/Biophilic/',imName{j},'.JPG'));

% Make image square 1536 x 1536 
targetSize = [1024 1024];
r = centerCropWindow2d(size(inImage), targetSize);
inImage=imcrop(inImage,r);

inImage = imresize(inImage, [667 667] )

imwrite(im2uint8(inImage),strcat('stimuli_color/cropsbio/',num2str(j),'.jpg'));

end
