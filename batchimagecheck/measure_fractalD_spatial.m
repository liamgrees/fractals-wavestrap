function [fractalD_spatial] = measure_fractalD_spatial(image_path,plot_flag)
%
% this code will measure the (spatial) fractal dimension (D) of your input
% image.
%
% indicate the path to your image (image_path) and whether or not you want
% plots to pop up (if yes, set plot_flag = 1).
%
% Log:
% 20210622 - Initialised. zji.

% the image is read and first its changed to grayscale image and then to
% black and white

if ischar((image_path)) ~= 1
    %     you probably directly loaded in an image...
    image = image_path;
else
    %     load image
    image=imread(image_path);
end

%% ensure input has equal dims
if mod(size(image,1),2) == 1 % will be 0 if even number, 1 if odd.
    
    %if size(f,1) is odd, then when it's divided by 2 later in the
    
    %script it will no longer be an integer... to ameliorate this, we
    
    %will scrape off one pixel in both x and y directions.
    
    image = image(1:end-1,:,:);
    
end

% y direction

if mod(size(image,2),2) == 1 % will be 0 if even number, 1 if odd.
    
    image = image(:,1:end-1,:);
    
end

% if size(image,1) ~= size(image,2)
%     
%     %     need to crop. will fix this in future versions
%     
%     s_original = [size(image,1) size(image,2)];
%     s_trim = min(pow2(fix(log2(s_original(1:2)))));
%     image = image(fix(s_original(1)/2) - s_trim/2+1:fix(s_original(1)/2) + s_trim/2, ...
%         fix(s_original(2)/2)-s_trim/2+1:fix(s_original(2)/2)+s_trim/2,:); %taken from one_over_f.
%     
% end

%%

orig_image = image;
if plot_flag == 1
    figure;subplot(1,3,1);imshow(orig_image);hold on;pause(0.25);title('input image')
end

image_size = size(image);

if numel(image_size) == 3; image=rgb2gray(image); end
orig_image=image;
image=image<mean(image(:));
[height,width,cols]=size(image);

if plot_flag == 1
    subplot(1,3,2);imshow(image);hold on;pause(0.25);title('thresholded image (@ mean luminance)')
end

% image = box_count_iterations_orig;
% initialises the varaibles for plotting the graph. scale is used to store
% the scaling factors and the count is used to store the number of boxes
% which contains parts of the image. For a given scale(1,i) value, the number
% of boxes occupied by the image will be available at count(1,i)

%% figure out how many scales can be used...

% scale=zeros(1,8);
% count=zeros(1,8);

% currently assuming imagehas even dimensions (ie is divisible by 2)...

%     scale = width/2;
%     sf = 2;

scale = gcd(height,width);
sf_width = width/scale;
sf_height = height/scale;

counter = 0;

while 1
    
    counter = counter + 1;
    
    if counter == 1
        
        newNum = scale/2;
        newScale_width = sf_width*2;
        newScale_height = sf_height*2;
        
    else
        
        newNum = newNum/2;
        newScale_width = newScale_width*2;
        newScale_height = newScale_height*2;
        
    end
    
    if floor(newNum)==newNum && newNum ~= 1
        
        scale = [scale newNum];
        sf_width = [sf_width newScale_width];
        sf_height = [sf_height newScale_height];
        
    else
        
        break
        
    end
    
end

if scale(1) == height && scale(1) == width
    scale = scale(2:end);
    sf_width = sf_width(2:end);
    sf_height = sf_height(2:end);
end

disp(['running boxcount across ' num2str(scale) ' scales'])

count = zeros(1,numel(scale));

% since the program is for a 512x512 image, the limit is set
% to 9 since 2^9=512
% if plot_orig_boxcount == 1;figure;end

box_plot = {};
grid_sizes_width = [];
grid_sizes_height = [];
piece_counter = [];

for i=1:numel(scale)
    
    %     parameters for box_plot variable: to be used for overlaying grid on
    %     original image for depiction purposes...
    
    box_plot{i} = zeros(height,width);
    
    % scaling factors are taken as 2,4,8,16... 512.
    
    % For each scaling factor, the total number of pieces are to be calculated,
    % and the number of pieces which contain the black dots (pixels) among them are to
    % be counted.
    
    % For eg, when the scaling factor is 2, it means the image is divided in to
    % half, hence we will get 4 pieces. And have to see how many of pieces
    % have the black dots.
    %     sf=2^i;
    %         pieces=sf(i)^2;
    
    pieceWidth=width/sf_width(i);
    pieceHeight=height/sf_height(i);
    pieces = (width/scale(i))*(height/scale(i));
    piece_counter = [piece_counter pieces];
    grid_sizes_width = [grid_sizes_width; pieceWidth];
    grid_sizes_height = [grid_sizes_height; pieceHeight];
    
    
    %initially we assume, we have 0 black pieces
    blackPieces=0;
    
    % Now we have to iterate through each pieces to see how many pieces have the
    % black dots (pixel) in it. We will consider the collection of pieces as
    % a matrix. We are counting from 0 for the ease of calculations.
    for pieceIndex=0:pieces-1
        
        % row and column indices of each pieces are calculated to estimate the
        % xy cordinates of the starting and ending of each piece.
        pieceRow=idivide(int32(pieceIndex), int32(sf_height(i))); %height
        pieceCol=rem(pieceIndex,sf_height(i));
        %             pieceCol=idivide(int32(pieceIndex), int32(sf_width(i))); %width
        
        xmin=(pieceRow*pieceHeight)+1;
        xmax=(xmin+pieceHeight)-1;
        ymin=(pieceCol*pieceWidth)+1;
        ymax=(ymin+pieceWidth)-1;
        
        % each piece is extracted and stored in another array for
        % convenience.
        
        eachPiece=image(ymin:ymax,xmin:xmax);
        %         each piece obtained is plotted on a plot for getting a view
        %         of the splitting of the whole image.
        
        if sum(eachPiece(:))~= size(eachPiece,1)*size(eachPiece,2) && ...
                sum(eachPiece(:))~= 0
            
            blackPieces=blackPieces+1;
            box_plot{i}(ymin:ymax,xmin:xmax) = ones(pieceHeight,pieceWidth);
            
        end
        
    end
    
    % the count of pieces which contains the black dots for a given scaling value
    % will be obtained here and will be stored in the respective variables.
    fprintf('%d\t->\t%d\n', scale(i), blackPieces);
    %     scale(1,i)=sf(i);
    count(1,i)=blackPieces;
end



% Now the process is over, the graph is plotted and the fractal dimension
% is calculated using the 'ployfit' function.

if range(count) == 0
    
    if plot_flag == 0; currFig = figure; hold on; else ; subplot(1,3,3); end
    
    disp('all boxes empty');
    
    fractalD_spatial = 0;
    disp(['reported fractal D: ' num2str(fractalD_spatial)])
    
    
%     plot(scale,count_for_log,'-o');
%     xlabel('scale')
%     ylabel('box count')
%     title('box count plot (linear axis)')
    
    if plot_flag == 0; close(currFig); end
    
else
    
    try
        
        if plot_flag == 0; currFig = figure; hold on; else ; subplot(1,3,3); end
        
        count_for_log = count;
        count_for_log(count_for_log == 0) = 1; %setting zeros to ones so logs work...
        
        plot(log(scale),log(count_for_log),'bo');
        params = polyfit(log(scale),log(count),1);
        y_params = polyval(params,log(scale));
        hold on;
        plot(log(scale),y_params,'r-');
        xlabel('scale')
        ylabel('box count')
        title('box count plot (linear axis)')
        fractalD_spatial = abs(params(1));
        disp(['reported fractal D: ' num2str(fractalD_spatial)])
        
        xlabel('log scale')
        ylabel('log box count')
        title(['fractal dim (linear fit): ' num2str(fractalD_spatial)]);
        
        if plot_flag == 0; close(currFig); end
        
    catch
        
        disp('not enough samples to determine fractal D... or something else went wrong.')
        close(gcf);
        
    end
    
end

%         figure;
%         plot(log(scale(1:end-1)),log(count(1:end-1)));
%         polyfit(log(scale(1:end-1)),log(count(1:end-1)),1);

%         figure;
%         plot(log(scale(2:end)),log(count(2:end)));
%        test_params = polyfit(log(scale(2:end)),log(count(2:end)),1);
%          disp(['test fractal D: ' num2str(test_params(1))])





%% plot overlays and draw grid lines
% test with an image with diff dims in x & y dimensions
if plot_flag == 1
    
    if numel(box_plot) > 6
        
        boxes_to_plot = 6;
        disp('only 6 box counts will be plotted to save computer memory');
        
    else
        
        boxes_to_plot = numel(box_plot);
        
    end
    
%     figure;
    
    for i = 1:boxes_to_plot
        %     for i = 1:(numel(box_plot)/2)
        
%         subplot(1,6,i);

        figure;
        
        h = gcf;
        
        %     B = labeloverlay(orig_image,box_plot{i},'Transparency',0.75);
        test = box_plot{i};
        B = labeloverlay(double(image),test,'Transparency',0.75);
        figure(h);imshow(B)
        
        hold on
        
        M = size(B,1);
        N = size(B,2);
        
        for k = 1:grid_sizes_width(i):M
            x = [1 N];
            y = [k k];
            %     plot(x,y,'Color','w','LineStyle','-');
            figure(h);plot(x,y,'Color','k','LineStyle','-','LineWidth',1);
        end
        
        for k = 1:grid_sizes_height(i):N
            x = [k k];
            y = [1 M];
            %     plot(x,y,'Color','w','LineStyle','-');
            figure(h);plot(x,y,'Color','k','LineStyle','-','LineWidth',1);
        end
        
        figure(h);title(['grid size: W -' num2str(grid_sizes_width(i)) ' H - ' num2str(grid_sizes_height(i)) ' pieces - ' num2str(count(i)) '/' num2str(piece_counter(i))]);
        
        hold off
        
        pause(0.75)
        
    end
    
end

end
