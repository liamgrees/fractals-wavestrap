%% createStimuli_color.m 
%  This script demonstrates how to apply wavestrapping to color images.
%  The procedure outlined here is described in section 2.3 of the following paper:
%  Puckett AM, Schira MM, Isherwood ZJ, Victor JD, Roberts JA, and 
%  Breakspear M. (2020) "Manipulating the structure of natural scenes using
%  wavelets to study the functional architecture of perceptual hierarchies 
%  in the brain" NeuroImage. 
%  https://www.sciencedirect.com/science/article/pii/S1053811920306595

%  Written by A.M.P. 2019.






N = 1536;

inImage = imread('stimuli_color/ZDB/test.JPG');

% Make image square 1536 x 1536 
targetSize = [1536 1536];
r = centerCropWindow2d(size(inImage), targetSize);
inImage=imcrop(inImage,r);

% Split out color components
inR_square = inImage(:,:,1);
inG_square = inImage(:,:,2);
inB_square = inImage(:,:,3);

%% Wavestrap each channel with the same seed --> will retain original color pallete across space

inR_wavC=wavesurr2(inR_square,[1 2 3 4 5 6 7 8],'db6',1,1);
inG_wavC=wavesurr2(inG_square,[1 2 3 4 5 6 7 8],'db6',1,1); 
inB_wavC=wavesurr2(inB_square,[1 2 3 4 5 6 7 8],'db6',1,1); 

%Adjust so histogram of wavestrapped iamge matches original
XR=inR_square(:);           
sX=inR_wavC(:);       %wavelet shuffled image
M(:,1)=sX; M(:,2)=1:N^2;
M=sortrows(M,1);
M(:,1)=sortrows(XR);
M=sortrows(M,2);
YRC=reshape(M(:,1),N,N);   %surrogate image with amplitude spectra of natural one

XG=inG_square(:);           
sX=inG_wavC(:);       %wavelet shuffled image
M(:,1)=sX; M(:,2)=1:N^2;
M=sortrows(M,1);
M(:,1)=sortrows(XG);
M=sortrows(M,2);
YGC=reshape(M(:,1),N,N);   %surrogate image with amplitude spectra of natural one

XB=inB_square(:);           
sX=inB_wavC(:);       %wavelet shuffled image
M(:,1)=sX; M(:,2)=1:N^2;
M=sortrows(M,1);
M(:,1)=sortrows(XB);
M=sortrows(M,2);
YBC=reshape(M(:,1),N,N);   %surrogate image with amplitude spectra of natural one

inR_wavC = uint8(YRC);
inG_wavC = uint8(YGC);
inB_wavC = uint8(YBC);

% Output wavestrapped image & crops
outImageC = cat(3,inR_wavC,inG_wavC,inB_wavC);

imwrite(inImage, 'stimuli_color/Crops/croptest.jpg')
imwrite(im2uint8(outImageC),'stimuli_color/Stimuliorder/testscramble.jpg');
