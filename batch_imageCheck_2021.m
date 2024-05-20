function compiledImageAnalysis = batch_imageCheck_2021(extension,plots,crop_type)
%
% This script will load all the images in your current folder and measure
% the following:
%
%   1) rms contrast
%   2) mean luminance
%   3) std luminance
%   4) fractal dimension
%   5) log-log fit of the 1/f amplitude spectrum (linear fit in log-log space)
%   6) linear fit of the 1/f amplitude spectrum (quadratic fit in linear space)
%
% This script needs to be run in the directory where all your
% images are. The extension of your images must be input as the first
% variable of this script, and you also have to indicate whether you want
% plots of the 1/f slope calculations.
%
% Examples of script input:
%
%   extension = 'bmp'
%   plots = 1 %input 0 for no plots
%   crop_type = 1
%                   1 - crop central region of image to largest poss. power of 2.
%                   2 - just ensure the dimensions are even. shave off 1 pixel in
%                       either direction to make image have even dims.
%                   3 - ensure dims are even, and crop by smallest dim in
%                       the centre to make the image square.
%
%   compiledImageAnalysis = batch_imageCheck('bmp',1,1)
%   compiledImageAnalysis = batch_imageCheck('tiff',0,2)
%
% Log:
% YYYYMMDD
% 20190410 - First created script. zji.
% 20210804 - Version 2. Updated 1/f slope and fractal D code.
%            Plots are currently not saved for the fractal D measurements.
% 20210805 - Version 3. Included a crop_type option where you can chose how
%            your images will be cropped. See above for more information 
%            of cropping methodologies.

timeScriptStarted = tic;
timeScriptStarted_forNaming = datestr(now, 30);

extension = 'jpg'
fractal_plots = 0;
crop_type = 2
plots = 0

% check crop type

if notDefined('crop_type')
    
    disp('crop_type not defined... defaulting to type 1: crop central region of image to largest possible power of 2');
    
    crop_type = 1;
    
elseif crop_type ~= 1 && crop_type ~= 2 && crop_type ~=3
    
    disp('crop_type not defined properly... defaulting to type 1: crop central region of image to largest possible power of 2');
    
    crop_type = 1;
    
end

% hide figures while script is running

set(0,'DefaultFigureVisible','off');

%% input variables

if ~isvar('extension')
    
    extension = ['bmp'];
    
end

if plots == 1
    
    spatialPlotsDir = [pwd '/' 'spatialPlots'];
    
    if ~isfolder(spatialPlotsDir)
        
        mkdir(spatialPlotsDir);
        
    end
    
end

%% load in images and get image size information

orig_image_size = {};
cropped_image_size = {};

images2process_dir = dir(['*' extension]);

% Get rid of ., .., and DS_Store

images2process_dir = images2process_dir(~ismember({images2process_dir.name}, ...
    {'.','..','.DS_Store','spatialPlots'}));

for ii = 1:numel(images2process_dir)
    
    images2process{ii} = imread(images2process_dir(ii).name);
    
    if size(images2process{ii},3) == 3
        
        images2process{ii} = rgb2gray(images2process{ii});
        
    end
    
    orig_image_size{ii} = size(images2process{ii});
    
    %% crop depending on crop_type
    
    if crop_type == 1 %crop central region of image to largest possible power of 2.
        
        if mod(size(images2process{ii},1),2) == 1 % will be 0 if even number, 1 if odd.
            
            %if size(f,1) is odd, then when it's divided by 2 later in the
            
            %script it will no longer be an integer... to ameliorate this, we
            
            %will scrape off one pixel in both x and y directions.
            
            images2process{ii} = images2process{ii}(1:end-1,:,:);
            
        end
        
        % y direction
        
        if mod(size(images2process{ii},2),2) == 1 % will be 0 if even number, 1 if odd.
            
            images2process{ii} = images2process{ii}(:,1:end-1,:);
            
        end
        
        if size(images2process{ii},1) ~= size(images2process{ii},2)
            
            %     need to crop. will fix this in future versions
            
            s_original = [size(images2process{ii},1) size(images2process{ii},2)];
            s_trim = min(pow2(fix(log2(s_original(1:2)))));
            images2process{ii} = images2process{ii}(fix(s_original(1)/2) - s_trim/2+1:fix(s_original(1)/2) + s_trim/2, ...
                fix(s_original(2)/2)-s_trim/2+1:fix(s_original(2)/2)+s_trim/2,:); %taken from one_over_f.
            
        end
        
    elseif crop_type == 2 %just ensure that dimensions are even
        
        if mod(size(images2process{ii},1),2) == 1 % will be 0 if even number, 1 if odd.
            
            %if size(f,1) is odd, then when it's divided by 2 later in the
            
            %script it will no longer be an integer... to ameliorate this, we
            
            %will scrape off one pixel in both x and y directions.
            
            images2process{ii} = images2process{ii}(1:end-1,:);
            
        end
        
        % y direction
        
        if mod(size(images2process{ii},2),2) == 1 % will be 0 if even number, 1 if odd.
            
            images2process{ii} = images2process{ii}(:,1:end-1);
            
        end
        
    elseif crop_type == 3 %ensure that the dimensions are even, and square (crop using smallest dim of image)
        
        
        if mod(size(images2process{ii},1),2) == 1 % will be 0 if even number, 1 if odd.
            
            %if size(f,1) is odd, then when it's divided by 2 later in the
            
            %script it will no longer be an integer... to ameliorate this, we
            
            %will scrape off one pixel in both x and y directions.
            
            images2process{ii} = images2process{ii}(1:end-1,:);
            
        end
        
        % y direction
        
        if mod(size(images2process{ii},2),2) == 1 % will be 0 if even number, 1 if odd.
            
            images2process{ii} = images2process{ii}(:,1:end-1);
            
        end
       
%         now make image square using the smallest dimension, and cropping
%         from the centre.
        
        [p3, p4] = size(images2process{ii});
        q1 = min([p3,p4]); %// size of the crop box
        i3_start = floor((p3-q1)/2); % or round instead of floor; using neither gives warning
        if i3_start == 0; i3_start = 1; end
        i3_stop = i3_start + q1 - 1;
        
        
      
        i4_start = floor((p4-q1)/2);
        if i4_start == 0; i4_start = 1; end
        i4_stop = i4_start + q1 - 1;
        
          
        
        
        images2process{ii} = images2process{ii}(i3_start:i3_stop,i4_start:i4_stop);
%         figure ,imshow(II);
        
        
    end
    
    cropped_image_size{ii} = size(images2process{ii});
    
end


%% image processing starts here...

imagesLoop_waitbar = waitbar(0,'Initializing waitbar...');

for imagesLoop = 1:numel(images2process)
    
    try
        
        waitbar((imagesLoop/numel(images2process)), ...
            imagesLoop_waitbar, ...
            ['Image ' num2str(imagesLoop) '/' ...
            num2str(numel(images2process)) ' being processed']);
        
    catch
        
        %         ignore...
        
    end
    
    % set up variables for storing info...
    
    rmsContrast = [];
    meanLum = [];
    stdLum = [];
    
    % S for spatial...
    
    S_fractalD = [];
    S_ampSlope = [];
    S_linear_slope = [];
    
    
    %% measure rms contrast
    
    [workingVar_rmsContrast,workingVar_meanLum,workingVar_stdLum] ...
        = calc_rmsContrast(images2process{imagesLoop});
    
    rmsContrast = workingVar_rmsContrast; %mean
    
    meanLum = workingVar_meanLum;
    
    stdLum = workingVar_stdLum; %std
    
    %% measure spatial fractal D
    
    %     workingVar = calc_fractalD_spatial(images2process{imagesLoop});
    
    workingVar = measure_fractalD_spatial(images2process{imagesLoop},fractal_plots);
    
    S_fractalD =  workingVar;
    
    %% measure spatial slope calculations
    
    if plots == 1
        
        [~,workingName,~] = fileparts(images2process_dir(imagesLoop).name);
        
        
        plotDir = spatialPlotsDir;
        
        plotName = [spatialPlotsDir '/' ...
            workingName '_spatialSlopePlot.png'];
        
    else
        
        plotDir = [];
        plotName = [];
        
    end
    
    %     [output_spatSlope] = calc_spatialSlope(images2process{imagesLoop},plots, ...
    %         plotDir,plotName);
    
    % [output_spatSlope] = calc_spatialSlope_R2(images2process{imagesLoop},plots, ...
    %     plotDir,plotName);

    
    % S_linear_slope = [output_spatSlope.spatialSlope_quadFit];
    % S_linear_R2 = [output_spatSlope.R2_quadFit];
    % 
    % S_ampSlope = [output_spatSlope.spatialSlope_logFit];
    % S_ampSlope_R2 = [output_spatSlope.R2_logFit];

    S_linear_slope = 0;
    S_linear_R2 = 0;
    
    S_ampSlope = 0;
    S_ampSlope_R2 = 0;
    
    %% save all variables into images2process variable...
    
    %imageName
    
    compiledImageAnalysis{imagesLoop}.fileName = images2process_dir(imagesLoop).name;
    
    % original image size
    
    compiledImageAnalysis{imagesLoop}.imageSize = size(images2process{imagesLoop});
    
    
    % cropped image size
    compiledImageAnalysis{imagesLoop}.croppedImageSize = size(images2process{imagesLoop});
    
    % rms contrast
    
    compiledImageAnalysis{imagesLoop}.rmsContrast = rmsContrast;
    
    % luminance
    
    compiledImageAnalysis{imagesLoop}.meanLum = meanLum;
    compiledImageAnalysis{imagesLoop}.stdLum = stdLum;
    
    % spatial fractal D
    
    compiledImageAnalysis{imagesLoop}.S_fractalD = S_fractalD;
    
    % log SS
    
    compiledImageAnalysis{imagesLoop}.logSS = S_ampSlope;
    compiledImageAnalysis{imagesLoop}.logSS_R2 = S_ampSlope_R2;
    
    % linear SS
    
    compiledImageAnalysis{imagesLoop}.linearSS = S_linear_slope;
    compiledImageAnalysis{imagesLoop}.linearSS_R2= S_linear_R2;
    
    %% save all variables into images2process variable...
    
    % save each loop...
    
    fileName_output = [timeScriptStarted_forNaming '_' 'raw' '_' ...
        'compiledImageAnalysis'];
    
    save(fileName_output,'compiledImageAnalysis');
    
    
end

close(imagesLoop_waitbar)

% transform output to something user friendly

%% compile together all data...
headers = {'Image Number', ... %1
    'Image Name', ...  %2
    'Original Image Size', ... %3
    'Cropped Image Size', ... %4
    'RMS Contrast', ...  %5
    'Spatial Mean Luminance', ...  %6
    'spatial STD Luminance', ...  %7
    'spatial fractal D', ...  %8
    'spatialSlope_logFit', ... %9
    'spatialSlope_logFit_R2', ... %10
    'spatialSlope_quadFit', ... %11
    'spatialSlope_quadFit_R2' }; %12

workingVar = cell((numel(compiledImageAnalysis)+1),numel(headers));

% headers

for ii = 1:numel(headers)
    
    workingVar{1,ii} = headers{ii};
    
end

for ii = 2:((numel(compiledImageAnalysis))+1)
    %image number 1
    
    workingVar{ii,1} = num2str(ii-1);
    
    % image names 2
    
    workingVar{ii,2} = compiledImageAnalysis{ii-1}.fileName;
    
    % original image size 3
    
    workingVar{ii,3} = orig_image_size{ii-1};%compiledImageAnalysis{ii-1}.imageSize;
    
    % cropped image size 4
    
    workingVar{ii,4} = cropped_image_size{ii-1}; %compiledImageAnalysis{ii-1}.imageSize;
    
    % contrast 5
    
    workingVar{ii,5} = compiledImageAnalysis{ii-1}.rmsContrast;
    
    % mean luminance 6
    
    workingVar{ii,6} = compiledImageAnalysis{ii-1}.meanLum;
    
    % std luminance 7
    
    workingVar{ii,7} = compiledImageAnalysis{ii-1}.stdLum;
    
    % fractal D 8
    
    workingVar{ii,8} = compiledImageAnalysis{ii-1}.S_fractalD;
    
    % 1/f slope - logFit 9
    
    workingVar{ii,9} = compiledImageAnalysis{ii-1}.logSS;
    
    % 1/f slope (r2) - logFit 10
    
    workingVar{ii,10} = compiledImageAnalysis{ii-1}.logSS_R2;
    
    % quad fit 11
    
    workingVar{ii,11} = compiledImageAnalysis{ii-1}.linearSS;
    
    % quad fit R2 12
    
    workingVar{ii,12} = compiledImageAnalysis{ii-1}.linearSS_R2;
    
end

compiledImageAnalysis = workingVar;

fileName_output = [timeScriptStarted_forNaming '_' 'userFriendly' '_' ...
    'compiledImageAnalysis' '_' 'cropType' '_' num2str(crop_type)];

save(fileName_output,'compiledImageAnalysis');

% remove working dir for extracted frames

set(0,'DefaultFigureVisible','on');

telapsed = toc(timeScriptStarted);
% minTime = min(telapsed,minTime);

disp(['This script took ' num2str(telapsed) ' minutes']);

end

function ndef = notDefined( varString )
    % Test whether a variable (usually a function argument) is defined
    %
    %    ndef = notDefined( varString )
    %
    % This routine is used to determine if a variable is defined in the calling
    % function's workspace.  A variable is defined if (a) it exists and (b) it
    % is not empty. This routine is used throughout the ISET code to test
    % whether arguments have been passed to the routine or a default should be
    % assigned.
    %
    % notDefined: 1 (true) if the variable is not defined in the calling workspace
    %             0 (false) if the variable is defined in the calling workspace
    %
    %  Defined means the variable exists and is not empty in the function that
    %  called this function.
    %
    %  This routine replaced many calls of the form
    %    if ~exist('varname','var') | isempty(xxx), ... end
    %
    %    with
    %
    %    if ieNotDefined('varname'), ... end
    %
    % bw summer 05 -- imported into mrVista 2.0
    % ras 10/05    -- changed variable names to avoid a recursion error.
    % ras 01/06    -- imported back into mrVista 1.0; why should we keep
    % typing 'ieNotDefined' in front of every function?

    if (~ischar(varString)), error('Varible name must be a string'); end

    ndef = 0;  % Assume the variable is defined

    str = sprintf('''%s''',varString);
    cmd1 = ['~exist(' str ',''var'') == 1'];
    cmd2 = ['isempty(',varString ') == 1'];
    cmd = [cmd1, ' | ',cmd2];

    % If either of these conditions holds, then not defined is true.
    ndef = evalin('caller',cmd1);     % Check that the variable exists in the caller space
    if ndef, return;                  % If it does not, return with a status of 0
    else
        ndef = evalin('caller',cmd2); % Check if the variable is empty in the caller space
        if ndef return;
        end
    end

end



