function [rmsContrast,meanLum,stdLum] = calc_rmsContrast(imageMatrix)

meanLum =  mean(imageMatrix(:));
stdLum = std(single(imageMatrix(:)));

rmsContrast = [stdLum/meanLum];

end
