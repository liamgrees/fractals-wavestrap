%% createStimuli_color.m 
%  This script demonstrates how to apply wavestrapping to color images.
%  The procedure outlined here is described in section 2.3 of the following paper:
%  Puckett AM, Schira MM, Isherwood ZJ, Victor JD, Roberts JA, and 
%  Breakspear M. (2020) "Manipulating the structure of natural scenes using
%  wavelets to study the functional architecture of perceptual hierarchies 
%  in the brain" NeuroImage. 
%  https://www.sciencedirect.com/science/article/pii/S1053811920306595

%  Written by A.M.P. 2019.


%% Input
imName = {'DSCN1572';...
    'DSCN2190';'DSCN2173';'DSCN2192';'DSCN1577';'DSCN2124';...
    'DSCN2123'; 'DSCN2172';'DSCN2184';'DSCN2156';'DSCN2078';'DSCN2179';...
    'DSCN2138';'DSCN2152';'DSCN2176';'DSCN2169'; 'DSCN1573'; 'DSCN1578';...
    'DSCN1579'; 'DSCN2087';'DSCN2101';'DSCN2141';'DSCN2174';'DSCN2144'; ...
    'DSCN1575';'DSCN1574';'DSCN1576';'DSCN2164'; 'DSCN2145'; 'DSCN2128';...
    'DSCN2119'; 'DSCN2189';'DSCN2122';'DSCN2126';'DSCN2135';'DSCN2137';...
    'DSCN2139';'DSCN2142';'DSCN2149';'DSCN2153'; 'DSCN2155'; 'DSCN2167'; ...
    'DSCN2170';'DSCN2175';'DSCN2181';'DSCN2183';'DSCN2185';'DSCN2187';...
    'DSCN2188';'DSCN2191'};

%% Input
for j = 1:49

N = 1536;

inImage = imread(strcat('stimuli_color/ZDB/',imName{j},'.JPG'));

% Make image square 1536 x 1536 
targetSize = [1536 1536];
r = centerCropWindow2d(size(inImage), targetSize);
inImage=imcrop(inImage,r);

% Split out color components
inR_square = inImage(:,:,1);
inG_square = inImage(:,:,2);
inB_square = inImage(:,:,3);

%% Wavestrap each channel with the same seed --> will retain original color pallete across space

inR_wavC=waveorder2(inR_square,[1 2 3 4 5 6 7 8 9 10],'db4',100,1);
inG_wavC=waveorder2(inG_square,[1 2 3 4 5 6 7 8 9 10],'db4',100,1); 
inB_wavC=waveorder2(inB_square,[1 2 3 4 5 6 7 8 9 10],'db4',100,1); 

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

imwrite(im2uint8(inImage),strcat('stimuli_color/Crops/crop_',num2str(j),'.jpg'));

imwrite(im2uint8(outImageC),strcat('stimuli_color/Stimuliorder/stimuli_',num2str(j),'.jpg'));

end
