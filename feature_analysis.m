close all;
clear;
clc;

% Browse file explorer to select an image
[filename, pathname] = uigetfile({'*.jpg;*.png;*.bmp','Image Files (*.jpg, *.png, *.bmp)'},'Select an Image');

% Check if user canceled file selection
if isequal(filename,0) || isequal(pathname,0)
    disp('User canceled file selection.');
    return;
end

% Load the selected image
fullpath = fullfile(pathname, filename);
k = imread(fullpath);

% Convert to grayscale
if size(k, 3) == 3
    k_gray = rgb2gray(k);
else
    k_gray = k;
end

% Call the function to compute and visualize edges
ked = computeAndVisualizeEdges(k_gray);

% Call the function to extract and visualize edge length
extractAndVisualizeEdgeLength(ked, k_gray);

% Count objects using connected component analysis
objectCount = countObjects(ked);

% Display the number of counted objects on the third figure
figure(3);
imshow(k);
hold on;
title(['Number of objects: ', num2str(objectCount)]);
hold off;

% Call the function to detect and extract corner features
detectAndVisualizeCorners(k_gray);

% Function to compute and visualize edges
function ked = computeAndVisualizeEdges(k_gray)
    % Convert to double for gradient computation
    k1 = double(k_gray);

    % Define Prewitt masks
    p_msk = [-1 0 1; -1 0 1; -1 0 1];

    % Compute gradients along x and y directions
    kx = conv2(k1, p_msk, 'same');
    ky = conv2(k1, p_msk', 'same');

    % Compute gradient magnitude
    ked = sqrt(kx.^2 + ky.^2);

    % Display the images
    figure;

    subplot(2, 2, 1);
    imshow(abs(kx), []);
    title('Edge Detection along X-axis');

    subplot(2, 2, 2);
    imshow(abs(ky), []);
    title('Edge Detection along Y-axis');

    subplot(2, 2, 3);
    imshow(abs(ked), []);
    title('Full Edge Detection');
end

% Function to extract and visualize edge length from the detected edges
function extractAndVisualizeEdgeLength(edgeImage, originalImage)
    % Extract edge segments
    [B, ~, N, ~] = bwboundaries(edgeImage);
    
    % Calculate edge length for each segment
    edgeLengths = zeros(N, 1);
    for k = 1:N
        boundary = B{k};
        edgeLengths(k) = size(boundary, 1);
    end
    
    % Display the original image with edge contours and annotations
    figure;
    imshow(originalImage);
    hold on;
    for k = 1:N
        boundary = B{k};
        plot(boundary(:, 2), boundary(:, 1), 'r', 'LineWidth', 2);
        text(mean(boundary(:, 2)), mean(boundary(:, 1)), num2str(edgeLengths(k)), 'Color', 'green', 'FontSize', 10);
    end
    hold off;
    title('Original Image with Edge Contours and Edge Length Annotations');
end

% Function to count objects using connected component analysis
function objectCount = countObjects(edgeImage)
    % Threshold the edge image to create a binary image
    binaryImage = edgeImage > 0.1; % You may need to adjust the threshold
    
    % Perform connected component analysis
    cc = bwconncomp(binaryImage);
    
    % Count the number of objects
    objectCount = cc.NumObjects;
end

% Function to detect and visualize corners
function detectAndVisualizeCorners(grayImage)
    % Detect Harris corners
    corners = detectHarrisFeatures(grayImage);
    
    % Extract features
    [features, valid_corners] = extractFeatures(grayImage, corners);
    
    % Display image
    figure;
    imshow(grayImage);
    hold on;
    
    % Plot valid corner points
    plot(valid_corners);
    
    hold off;
    title('Image with Harris Corner Features');
end
