function [output] = calc_spatialSlope_R2(input,plotFlag,plotDir,plotName,save_output)

% principle code written by cc.
% batch script put together by zji.
%
% Data spit out from this function:
%
%   1) output.spatialSlope_quadFit
%   2) output.spatialSlope_logFit;
%   3) output.R2_quadFit;
%   4) output.R2_logFit;
%   5) output.chosenCutOff
%   6) output.maxPointsTested
%   7) output.data_and_fit_loglog = [];
%   8) output.data_and_fit_linear_linearaxis
%
% Log:
%
% YYYYMMDD
% 20181031 - Version 2 created.
% 20210504 - Version 3. Added data + fit output in matrix form for plotting
% outside of matlab (7 to 9).
%
% Includes new method of fitting.. This script will exclude a certain
% number of points at the beginning of the data set to account for low SFs
% biasing fits. The number of points chosen to be removed is based on which
% number has the lowest AIC.
%
% This script will also allow non-square input... However, non-square (and
% also input that isn't a factor of 2^n) takes longer to process since
% Fourier transforms will take longer.
%
% This script will also check if your input is coloured and covert to
% grayscale if it is.
%
% To do for future versions:
%
% -Clean up plotting section of the code. A lot of redundancies.
%
% zji.

%% check if we should run this in parallel

if ~isvar('parallel_flag')

    %save in the current directory in a directory called plots

    parallel_flag = 1;


elseif isempty(parallel_flag)

    %save in the current directory in a directory called plots

    parallel_flag = 1;


end
%
% if parallel_flag == 1
%
%     parforArg = Inf; %use all cores
%     disp('running in paralell...')
% else
%     parforArg = 0; %use only 1 core.
%     disp('running in serial...')
% end

%% input vars

scriptStartTime = datestr(now, 30);

cutOffPercentage = 0.02; %try fitting w/ and w/o bottom 2% of points
% sometimes there are large magnitude biases across low frequencies which
% can greatly affect the fit of the spectrum in linear space. As such, we
% try to fit w/ and w/o the bottom indicated percentage of points to see if
% we can achieve a better fit. Originally this value was the first '10'
% points, but this wouldn't work for all dataset... so we've opted for a
% percentage going forwared.


%% input checks

% if ~isfolder(movie_path{1})
%
%     [input] = double(load_frames_from_movie(movie_path{1}));
%
% else % you've input a folder full of frames...?
%
%     [input] = load_frames(movie_path{1},movie_path{2},1);
%
% end


frames = size(input,3);

% maxCutOff = ceil(frames*cutOffPercentage); %cut of bottom 2% of points


% will the input be an integer when divided by 2? if not, remove a pixel in
% the x and y direction.

% x direction

if mod(size(input,1),2) == 1 % will be 0 if even number, 1 if odd.

    %if size(f,1) is odd, then when it's divided by 2 later in the

    %script it will no longer be an integer... to ameliorate this, we

    %will scrape off one pixel in both x and y directions.

    input = input(1:end-1,:,:);

end

% y direction

if mod(size(input,2),2) == 1 % will be 0 if even number, 1 if odd.

    input = input(:,1:end-1,:);

end

if size(input,1) ~= size(input,2)
    %
    %     %     need to crop. will fix this in future versions
    %
    %     s_original = [size(input,1) size(input,2)];
    %     s_trim = min(pow2(fix(log2(s_original(1:2)))));
    %     input = input(fix(s_original(1)/2) - s_trim/2+1:fix(s_original(1)/2) + s_trim/2, ...
    %         fix(s_original(2)/2)-s_trim/2+1:fix(s_original(2)/2)+s_trim/2,:); %taken from one_over_f.

    % make image square by padding:

    [input,~] = makeimagesquare(input);

end

maxCutOff = ceil(size(input,1)*cutOffPercentage);

%% calculate spatial slope across frames

spatialSlope_quadFit = zeros(1,frames);
spatialSlope_logFit = zeros(1,frames);
R2_quadratic = zeros(1,frames);
R2_loglog = zeros(1,frames);
chosenCutOff = zeros(1,frames);
A_loglog = {};
B_loglog = {};
f_loglog = {};
A_quadratic = {};
B_quadratic = {};
f_quadratic = {};

data_and_fit_linear_linearaxis = {};
data_and_fit_loglog = {};

% parfor (ii = 1:frames,parforArg)
for ii = 1:frames
    frame = double(input(:,:,ii));
    % subtract DC component
    m = mean(frame(:));
    frame = frame - m;

    %     frame = imcomplement(frame);

    % compute its 2-D spectrum

    spec = fftn(frame);

    % convert from complex to amplitude spectrum

    ampspec = abs(spec);

    % collapse across spatial dimensions

    x = mean(ampspec,1);
    xx = mean(x,2);

    xsize = size(frame,1);
    ysize = size(frame,2);
    imf=fftshift(fft2(frame));

    impf=abs(imf);

    Pf = rotavg(impf);

    x_vec = 1:xsize/2;
    y_vec = Pf(2:(ysize/2)+1);
    y_vec = y_vec';

    %% figure out what cutoff gives the best fit

    aic_quadratic = [];
    aic_loglog = [];

    r2_comparison_quad = [];
    r2_comparison_loglog = [];


    %%

    for cutOffTest_loop = 1:maxCutOff

        % linear

        % fit hyperbolic to raw data ...
        maxParam = max(y_vec(cutOffTest_loop:end));
        init_params = [maxParam -1]; % first guess to start iteration ...

        %         try
        params_0 = nlinfit(x_vec(cutOffTest_loop:end),y_vec(cutOffTest_loop:end), ...
            'one_over_x',init_params);

        %         catch
        %
        %             %perhaps the first value is messing things up
        %
        %             maxParam = max(y_vec(3:end));
        %             init_params = [maxParam -1]; % first guess to start iteration ...
        %
        %
        %             params_0 = nlinfit(x_vec(3:end),y_vec(3:end), ...
        %                 'one_over_x',init_params);
        %
        %
        %         end

        spatialSlope_quadFit(ii) = params_0(2);

        probs = one_over_x(params_0,x_vec);

        % loglog

        A_loglog{ii} = log(x_vec(cutOffTest_loop:end));
        B_loglog{ii} = log(y_vec(cutOffTest_loop:end));
        b_loglog = polyfit(A_loglog{ii}, B_loglog{ii}, 1);

        spatialSlope_logFit(ii) = b_loglog(1);


        %% Calculate R2

        % quad fit

        A_quadratic{ii} = x_vec(cutOffTest_loop:end);
        B_quadratic{ii} = y_vec(cutOffTest_loop:end);
        f_quadratic{ii} = probs(cutOffTest_loop:end);
        Bbar_quadratic = mean(B_quadratic{ii});
        SStot_quadratic = sum((B_quadratic{ii} - Bbar_quadratic).^2);
        SSres_quadratic = sum((B_quadratic{ii} - f_quadratic{ii}).^2);
        R2_quadratic(ii) = 1 - SSres_quadratic/SStot_quadratic;

        quadratic_residuals = B_quadratic{ii} - f_quadratic{ii};

        % linear fit (logged)

        f_linear = polyval(b_loglog, A_loglog{ii});
        f_linear = exp(f_linear);
        Bbar_linear = mean(B_quadratic{ii});
        SStot_linear = sum((B_quadratic{ii} - Bbar_linear).^2);
        SSres_linear = sum((B_quadratic{ii} - f_linear).^2);
        R2_linearSpace = 1 - SSres_linear/SStot_linear;

        linear_residuals = B_quadratic{ii} - f_linear;


        %         incorrect placement of the following.. we want to compare the log
        %         log fit to the quad fit!!!


        %         workingVar_aic_quadratic = aic(3, numel(x_vec), quadratic_residuals);
        %         workingVar_aic_linear = aic(2, numel(x_vec),linear_residuals);
        %
        %         aic_quadratic = [aic_quadratic; workingVar_aic_quadratic ...
        %             workingVar_aic_linear];
        %
        %         r2_comparison_quad = [r2_comparison_quad; R2_quadratic(ii) R2_linearSpace];



        %% loglog fit

        A_loglog{ii} = log(x_vec);
        B_loglog{ii} = log(y_vec);
        b_loglog = polyfit(A_loglog{ii}, B_loglog{ii}, 1);

        f_loglog{ii} = polyval(b_loglog, A_loglog{ii});
        Bbar_loglog = mean(B_loglog{ii});
        SStot_loglog = sum((B_loglog{ii} - Bbar_loglog).^2);
        SSreg_loglog = sum((f_loglog{ii} - Bbar_loglog).^2);
        SSres_loglog = sum((B_loglog{ii} - f_loglog{ii}).^2);
        R2_loglog(ii) = 1 - SSres_loglog/SStot_loglog;
        loglog_residuals = f_loglog{ii} - Bbar_loglog;

        % aic_loglog = [aic_loglog; aic([2, numel(x_vec), loglog_residuals)];
         aic_loglog = [aic_loglog; aicbic(loglog_residuals,2, numel(x_vec))];
        r2_comparison_loglog = [r2_comparison_loglog; R2_loglog(ii)];

        %         compare loglog to quad

        % workingVar_aic_quadratic = aic(2, numel(x_vec), quadratic_residuals);
        % workingVar_aic_linear = aic(1, numel(x_vec),loglog_residuals);

        workingVar_aic_quadratic = aicbic(quadratic_residuals,2,numel(x_vec));
        workingVar_aic_linear = aicbic(loglog_residuals,1,numel(x_vec));

        aic_quadratic = [aic_quadratic; workingVar_aic_quadratic' ...
            workingVar_aic_linear'];

        r2_comparison_quad = [r2_comparison_quad; R2_quadratic(ii) R2_loglog(ii)];

        %probs with loglog slope...

        %         probs_loglog = one_over_x([params_0(1) ],x_vec);

        if cutOffTest_loop == maxCutOff

            %% get point before inflection

            %             summaryDiffs_AIC = (abs(aic_quadratic(:,1)-aic_quadratic(:,2)));
            %             summaryDiffs_AIC = (abs(aic_quadratic(:,2)-aic_quadratic(:,1)));

            %             chosenCutOff(ii) = find(summaryDiffs_AIC == min(summaryDiffs_AIC))-1;
            
            try
            chosenCutOff(ii) = find(r2_comparison_quad(:,1) == max(r2_comparison_quad(:,1))); %find highest R2 fit

            if chosenCutOff(ii) == 0

                %this value means that you're better off not getting rid of any
                %points.

                chosenCutOff(ii) = 1;

            end

            catch

                disp('r2 comparison failed... not cutting off any points!')

                chosenCutOff(ii) = 1;

            end

            % recalc the slope and other variables for plotting for the
            % determined cutOff

            % linear

            % fit hyperbolic to raw data ...
            maxParam = max(y_vec(chosenCutOff(ii):end));
            init_params = [maxParam -1]; % first guess to start iteration ...


            %             try
            params_0 = nlinfit(x_vec(chosenCutOff(ii):end),y_vec(chosenCutOff(ii):end), ...
                'one_over_x',init_params);

            %             catch
            %
            %                 %                 maybe the first value is throwing things off.
            %
            %                 maxParam = max(y_vec(2:end));
            %                 init_params = [maxParam -1]; %
            %                 params_0 = nlinfit(x_vec(2:end),y_vec(2:end), ...
            %                     'one_over_x',init_params);
            %
            %
            %             end

            spatialSlope_quadFit(ii) = params_0(2);

            probs = one_over_x(params_0,x_vec);

            % loglog

            A_loglog{ii} = log(x_vec(chosenCutOff(ii):end));
            B_loglog{ii} = log(y_vec(chosenCutOff(ii):end));
            b_loglog = polyfit(A_loglog{ii}, B_loglog{ii}, 1);

            spatialSlope_logFit(ii) = b_loglog(1);


            %% Calculate R2

            % quad fit

            A_quadratic{ii} = x_vec(chosenCutOff(ii):end);
            B_quadratic{ii} = y_vec(chosenCutOff(ii):end);
            f_quadratic{ii} = probs(chosenCutOff(ii):end);
            Bbar_quadratic = mean(B_quadratic{ii});
            SStot_quadratic = sum((B_quadratic{ii} - Bbar_quadratic).^2);
            SSres_quadratic = sum((B_quadratic{ii} - f_quadratic{ii}).^2);
            R2_quadratic(ii) = 1 - SSres_quadratic/SStot_quadratic;

            data_and_fit_linear_linearaxis{ii} = [A_quadratic{ii}'  B_quadratic{ii}' f_quadratic{ii}'];



            % linear fit (data fit on loglog axis, transformed to linear axis)

            f_linear = polyval(b_loglog, A_loglog{ii});
            f_linear = exp(f_linear);
            Bbar_linear = mean(B_quadratic{ii});
            SStot_linear = sum((B_quadratic{ii} - Bbar_linear).^2);
            SSres_linear = sum((B_quadratic{ii} - f_linear).^2);
            R2_linearSpace = 1 - SSres_linear/SStot_linear;

            %% loglog fit

            A_loglog{ii} = log(x_vec(chosenCutOff(ii):end));
            B_loglog{ii} = log(y_vec(chosenCutOff(ii):end));
            b_loglog = polyfit(A_loglog{ii}, B_loglog{ii}, 1);

            spatialSlope_logFit(ii) = b_loglog(1);

            f_loglog{ii} = polyval(b_loglog, A_loglog{ii});
            Bbar_loglog = mean(B_loglog{ii});
            SStot_loglog = sum((B_loglog{ii} - Bbar_loglog).^2);

            SSres_loglog = sum((B_loglog{ii} - f_loglog{ii}).^2);
            R2_loglog(ii) = 1 - SSres_loglog/SStot_loglog;

            data_and_fit_loglog{ii} = [A_loglog{ii}' B_loglog{ii}' f_loglog{ii}'];

        end
    end % cutOffTest_loop = 1:maxCutOff


    if save_output == 1
        if ~isvar('plotDir')

            %save in the current directory in a directory called plots

            dirForSaving = ['spatSlopeAnalysis'];


            if ~isdir(dirForSaving)

                mkdir(dirForSaving)

            end

        elseif isempty(plotDir)

            %save in the current directory in a directory called plots

            dirForSaving = ['spatSlopeAnalysis'];

            mkdir(dirForSaving)

            if ~isdir(dirForSaving)

                mkdir(dirForSaving)

            end

        else

            dirForSaving = plotDir;

            if ~isdir(dirForSaving)

                mkdir(dirForSaving)

            end

        end

    end


    if plotFlag == 1

        % loglog fit

        figure;

        subplot 131

        plot(A_loglog{ii},B_loglog{ii},'bo');

        hold on

        plot([A_loglog{ii}(1), A_loglog{ii}(end)], b_loglog(2) ...
            + b_loglog(1).*[A_loglog{ii}(1),A_loglog{ii}(end)],'r-');

        hold off

        title(['log fit on log axis']);
        xlabel('log SF');
        ylabel('log amp');

        legend(['slope ' num2str(spatialSlope_logFit(ii))],['R^2 ' ...
            num2str(R2_loglog(ii))]);

        % linear fit on linear axis

        subplot 132

        hold on;

        plot(A_quadratic{ii},B_quadratic{ii},'bo'); % data
        hold on;

        plot(A_quadratic{ii},f_quadratic{ii},'r-'); % fit of a*(x^b)
        hold off;

        title(['linear fit on linear axis']);
        xlabel('SF');
        ylabel( 'amp');

        legend(['slope ' num2str(spatialSlope_quadFit(ii))],['R^2 ' num2str(R2_quadratic(ii))]);

        % linear fit on loglog axis

        A_linearFitLogAxis = log(x_vec(chosenCutOff(ii):end));
        B_linearFitLogAxis = log(y_vec(chosenCutOff(ii):end));
        f_linearFitLogAxis = log(probs(chosenCutOff(ii):end));
        Bbar_linearFitLogAxis = mean(B_linearFitLogAxis);
        SStot_linearFitLogAxis = sum((B_linearFitLogAxis - Bbar_linearFitLogAxis).^2);
        SSres_linearFitLogAxis = sum((B_linearFitLogAxis - f_linearFitLogAxis).^2);
        R2_linearFitLogAxis = 1 - SSres_linearFitLogAxis/SStot_linearFitLogAxis;

        subplot 133

        plot(A_linearFitLogAxis,B_linearFitLogAxis,'bo');

        hold on;

        plot(A_linearFitLogAxis,f_linearFitLogAxis,'r-'); % fit of a*(x^b)

        hold off;

        title(['linear fit on logged axis']);
        xlabel('log SF');
        ylabel('log amp');

        legend(['slope ' num2str(spatialSlope_quadFit(ii))],['R^2 ' num2str(R2_linearFitLogAxis)]);

        % save figures



        if ~isvar('plotName')

            plotName = [dirForSaving '/' 'spatialSlopePlot_' scriptStartTime ...
                '_cutOff' num2str(chosenCutOff(ii)) '.png'];

        elseif  isempty(plotName)

            plotName = [dirForSaving '/' 'spatialSlopePlot_' scriptStartTime ...
                '_cutOff' num2str(chosenCutOff(ii)) '.png'];

        end

        if save_output == 1
            saveas(gcf,plotName)
        end
        %         close(gcf)

    end

    % prepare everything for output...



end

output = struct;

output.spatialSlope_quadFit = spatialSlope_quadFit;
output.spatialSlope_logFit = spatialSlope_logFit;
output.R2_quadFit = R2_quadratic;
output.R2_logFit = R2_loglog;
output.chosenCutOff = chosenCutOff;
output.maxPointsTested = maxCutOff;
output.data_and_fit_loglog = data_and_fit_loglog; %[A_loglog' B_loglog' f_loglog'];
output.data_and_fit_linear_linearaxis = data_and_fit_linear_linearaxis; %[A_quadratic B_quadratic f_quadratic];
% output.data_and_fit_linear_loglogaxis = [log(A_quadratic)' log(B_quadratic)' log(f_quadratic)'];
% save output

if save_output == 1

    dataSaveName = [dirForSaving '/' 'spatialSlopeOutput_' scriptStartTime '.mat'];

    save(dataSaveName,'output')

end

end

function probs = one_over_x(params,x)

% Design matrix for f(x) = a*x^b.
%
%
% % CC 6.3.01

offSetParam = 0; %100; %fit rounding later down the track...

probs = params(1).*(x.^(params(2)))+offSetParam;

end

function [outputim,maskmap] = makeimagesquare(input)
% pads your image with zeros if it's not square
% eg, a 390x512 image would become 512x512
% if your image is already square, the code will break.
% image needs to have even dimensions and also square (x = y)
% this code will also spit out a mask map so the image can be easily
% cropped after analysis.
%
% can use maskmap to crop image back to original size with following lines;
%
%   [row_start cols_start] = find(maskmap, 1, 'first')
%   [row_end cols_end] = find(maskmap, 1, 'last')
%   finalim = outputim(row_start:row_end,cols_start:cols_end);
%
% Log:
% 20200311. Initialised. zoeyisherwood.
% Contact: zoey.isherwood@gmail.com

% check dims---------------------------------------------------------------

if ndims(input) == 3
    disp(['your input is colored. your input will be' ...
        ' converted to greyscale...'])
    input = input(:,:,1:3); %just in case there's an alpha channel in there.
    input = rgb2gray(input);

end

% get sizes----------------------------------------------------------------

a = size(input,1);
b = size(input,2);

if a == b %end function if the image is already square

    disp('your input is already square.')
    outputim = input;
    maskmap = ones(a,b);

else

    disp('your image is not square. converting...');

    maskmap = ones(a,b);

    %determine padding amount----------------------------------------------

    %firstly, which dim is larger?

    padsize = max(a,b);

    % ... and is this dim an even number? if not, add 1.

    if mod(padsize,2) ~= 0 %it's an odd number

        disp(['large dim isn''t even. padding image to make dims even'])

        padsize = padsize + 1;

    end

    % pad depending on dimension-------------------------------------------

    outputim = zeros(padsize,padsize);

    start_row = ceil((padsize-a)/2 + 1);
    start_col = ceil((padsize-b)/2 + 1);

    outputim(start_row:start_row+a-1,start_col:start_col+b-1) = input;

    maskmap = zeros(padsize,padsize);
    maskmap(start_row:start_row+a-1,start_col:start_col+b-1) = ones(a,b);

    disp('done')

end

end
