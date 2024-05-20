%% Input
imName = {'1';...
    '2';'3';'4';'5';'6';...
    '7'; '8';'9';'10';'11';'12';...
    '13';'14';'15'};

%% Input
for j = 1:15


inImage = imread(strcat('stimuli_color/cropsbiophilic/',imName{j},'.jpg'));

% Make image square 1536 x 1536 
inImage = imresize(inImage, [667 667] )

imwrite(im2uint8(inImage),strcat('stimuli_color/cropstransformedbio/',num2str(j),'.jpg'));

end