  
    %Read image in, convert to grayscale, and make square. 
    inImage = imread('DSCN1572.JPG');
    inImage = rgb2gray(inImage);
    imageOrig = im2double(inImage);

       
    % B2_2
    imageN1S1F1R1 = waveorder2(imageOrig,[1 2 3 4 5 6 7 8 9 10], 'db4', 100,1); % wavestrap using rectsurr2
    % rectsur2 resampled a central rectangular region. Here we are resampling the
    % inner 19/20ths while leaving a 1/20 border untouched. Then we crop.
    % This approach was chosen due to issues with edge artifacts. 
   
    imageN1S1F1R1Resize = imageN1S1F1R1(86:1100,86:1100); %Crop
    imageN1S1F1R1Resize = imresize(imageN1S1F1R1Resize,0.5835); %Now resize 

    %Write out image and clear some variables
    imwrite(imageN1S1F1R1Resize,'wavestrapped/DSCN1572.jpg')
    clear sX M Y

