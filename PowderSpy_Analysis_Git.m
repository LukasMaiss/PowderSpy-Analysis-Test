clc, clear, close all

%% Parameters for Analysis - DO NOT CHANGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
pixel_size = 0.02327;
x_offset = 500;
y_offset = 1300;
nozzle_offset = 50;
threshold = 10;
cog_threshold = 100;
intensity_z_threshold = 1100;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Parameters for evaluation analysis
focus_1D_86 = 1;
focus_2D_86 = 1;
focus_fwhm = 1;

focus_top = 850;
focus_bottom = 750;
%% ONLY FOR MULTIJET NOZZLES
multijet_nozzle = 0; % if nozzle to analyse is multijet nozzle set value to 1
number_injectors = 6;
threshold_multijet = 25;

%% load image stack
processedImages = readNPY('processedImages.npy');

%% removing single particles
for i=1:size(processedImages,3)
    bw = processedImages(:,:,i);
    bw(bw<10)=0;
    bw(bw>0)=1;
    bw_cut=bwareaopen(bw, 30, 4);
    copy = processedImages(:,:,i);
    copy(bw_cut==0)=0;
    processedImages(:,:,i)=copy;
end

%% visualize image stack
% figure
for i=1:10:1300
%     imshow(processedImages(:,:,1301-i), Colormap=jet)
%     disp(i)
end

%% underlaying matrices for report visualization
xz_cross = uint8(zeros(size(processedImages,2),size(processedImages,3),size(processedImages,1)));

for i=1:size(xz_cross,3)
    xz_cross(:,:,i) = processedImages(i,:,:);
end

yz_cross = uint8(zeros(size(processedImages,1),size(processedImages,3),size(processedImages,2)));
for i=1:size(yz_cross,3)
    yz_cross(:,:,i) = processedImages(:,i,:);
end

%% working distance and powder focus 1D 86% rule

if focus_1D_86 == 1
    [x,y,z] = ind2sub(size(processedImages), find(processedImages>240));
    
    z_ebene = round((min(z)+max(z(z<intensity_z_threshold)))/2);

    check_focus_layer = isempty(z_ebene);

    if check_focus_layer == 0
        focus_start = z_ebene;
        focus_end = z_ebene;

        diameter1_1D_86 = zeros(41,1);
        ratio1_1D_86 = zeros(41,1);
        diameter2_1D_86 = zeros(41,1);
        ratio2_1D_86 = zeros(41,1);
        diameter_1D_86 = zeros(1,1);
    else
        focus_start = focus_bottom;
        focus_end = focus_top;

        diameter1_1D_86 = zeros(41,focus_end - focus_start);
        ratio1_1D_86 = zeros(41,focus_end - focus_start);
        diameter2_1D_86 = zeros(41,focus_end - focus_start);
        ratio2_1D_86 = zeros(41,focus_end - focus_start);
        diameter_1D_86 = zeros(focus_end - focus_start,1);
    end

    l=0;
    for k=focus_start:focus_end

        layer = processedImages(:,:,k);
        layer(layer<threshold) = 0;
        cog_layer = layer;
        cog_layer(cog_layer<cog_threshold) = 0;

        binaryImage = true(size(cog_layer));
        measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
        cog = measurements(1).WeightedCentroid;
        max_y = round(cog(1));
        max_x = round(cog(2));

        j=0;
        l=l+1;
        for i=max_x-20:max_x+20
            profil1 = double(layer(i,:));
            int_profil1 = sum(profil1);
            int_86_1 = 0;
            j=j+1;
            rad1=1;
            max_profil1 = max_y;

            while int_86_1<0.86*int_profil1
                int_86_1=sum(profil1(max_profil1-rad1:max_profil1+rad1));
                rad1=rad1+1;
            end

            diameter1_1D_86(j,l) = 2*rad1 + 1;
            ratio1_1D_86(j,l) = diameter1_1D_86(j,l)/max(profil1);
        end
        j=0;
        for i=floor(mean(max_y))-20:floor(mean(max_y))+20
            profil2 = double(layer(:,i));
            int_profil2 = sum(profil2);
            int_86_2 = 0;
            j=j+1;
            rad2=1;
            max_profil2 = max_x;

            while int_86_2<0.86*int_profil2
                int_86_2=sum(profil2(max_profil2-rad2:max_profil2+rad2));
                rad2=rad2+1;
            end

            diameter2_1D_86(j,l) = 2*rad2 + 1;
            ratio2_1D_86(j,l) = diameter2_1D_86(j,l)/max(profil2);
        end

        [x_V1, y_V1] = find(ratio1_1D_86(:,l) == min(ratio1_1D_86(:,l)));
        min_durchmesser_x = diameter1_1D_86(x_V1, l);
        [x_V2, y_V2] = find(ratio2_1D_86(:,l) == min(ratio2_1D_86(:,l)));
        min_durchmesser_y = diameter2_1D_86(x_V2, l);
        diameter_1D_86(l) = (min_durchmesser_y(1) + min_durchmesser_x(1))/2;
    end
    [val_1D_86,focusz_val] = min(diameter_1D_86);
    focusz_pos = find(diameter_1D_86==val_1D_86);
    focusz = max(focusz_pos);
    fd_1D_86 = round(val_1D_86*pixel_size,2);
    wd_1D_86 = round((y_offset-focus_start-focusz+nozzle_offset)*pixel_size,2);

    disp(['Focus diameter 1D 86%: ',num2str(fd_1D_86), ' mm; Working distance 1D 86%: ', num2str(wd_1D_86), ' mm' ])

    focus_position = y_offset - (y_offset - focus_start - focusz);
    focus_layer = processedImages(:,:,focus_position);

end

cog_layer = focus_layer;
cog_layer(cog_layer<cog_threshold) = 0;
binaryImage = true(size(cog_layer));
measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
cog = measurements(1).WeightedCentroid;
max_y = round(cog(1));
max_x = round(cog(2));
figure
imshow(focus_layer, Colormap=jet)
viscircles([max_y max_x], round(val_1D_86/2));

%% working distance and powder focus 2D 86% rule

if focus_2D_86 == 1
    [x,y,z] = ind2sub(size(processedImages), find(processedImages>240));
    
    z_ebene = round((min(z)+max(z(z<intensity_z_threshold)))/2);

    check_focus_layer = isempty(z_ebene);

    if check_focus_layer == 0
        focus_start = z_ebene;
        focus_end = z_ebene;

        results_2D_86 = zeros(focus_end-focus_start,3);
    else
        focus_start = focus_bottom;
        focus_end = focus_top;

        results_2D_86 = zeros(focus_end-focus_start,3);
    end

    l=0;
    for k=focus_start:focus_end

        layer = processedImages(:,:,k);
        layer(layer<threshold) = 0;
        cog_layer = layer;
        cog_layer(cog_layer<cog_threshold) = 0;

        binaryImage = true(size(cog_layer));
        measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
        cog = measurements(1).WeightedCentroid;
        max_y = round(cog(1));
        max_x = round(cog(2));


        sum_intensity = layer(max_x,max_y);
        j = 0;

        [columnsInImage, rowsInImage] = meshgrid(1:size(layer,2), 1:size(layer,1));
        while sum_intensity <= 0.86 * sum(layer, 'all')
            j=j+1;
            radius = j;
            circlePixels = (rowsInImage - max_x).^2 + (columnsInImage - max_y).^2 <= radius.^2;
            sum_intensity = sum(double(layer) .* double(circlePixels),'all');

        end
        focusdiameter = 2*j+1;
        maxintensity = max(layer,[],[1 2]);

        results_2D_86(k-focus_start+1,1) = double(focusdiameter) / double(maxintensity);
        results_2D_86(k-focus_start+1,2) = focusdiameter;
        results_2D_86(k-focus_start+1,3) = maxintensity;
    end
    [val,focusz] = min(results_2D_86(:,1));
    fd_2D_86 = round(results_2D_86(focusz, 2)*pixel_size,2);
    wd_2D_86 = round((y_offset-focus_start-focusz+nozzle_offset)*pixel_size,2);

    disp(['Focus diameter 2D 86%: ',num2str(fd_2D_86), ' mm; Working distance 2D 86%: ', num2str(wd_2D_86), ' mm' ])

    focus_position = y_offset - (y_offset - focus_start - focusz);
    focus_layer = processedImages(:,:,focus_position);

end

cog_layer = focus_layer;
cog_layer(cog_layer<cog_threshold) = 0;
binaryImage = true(size(cog_layer));
measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
cog = measurements(1).WeightedCentroid;
max_y = round(cog(1));
max_x = round(cog(2));
figure
imshow(focus_layer, Colormap=jet)
viscircles([max_y max_x], round(results_2D_86(focusz,2)/2));

%% working distance and powder focus 1D full width at half maximum

if focus_fwhm == 1
    [x,y,z] = ind2sub(size(processedImages), find(processedImages>240));

    z_ebene = round((min(z)+max(z(z<intensity_z_threshold)))/2);

    check_focus_layer = isempty(z_ebene);

    if check_focus_layer == 0
        focus_start = z_ebene;
        focus_end = z_ebene;

        diameter1_fwhm = zeros(41,1);
        ratio1_fwhm = zeros(41,1);
        diameter2_fwhm = zeros(41,1);
        ratio2_fwhm = zeros(41,1);
        ratio_fwhm = zeros(1,1);
        diameter_fwhm = zeros(1,1);
    else
        focus_start = focus_bottom;
        focus_end = focus_top;

        diameter1_fwhm = zeros(41,focus_end - focus_start);
        ratio1_fwhm = zeros(41,focus_end - focus_start);
        diameter2_fwhm = zeros(41,focus_end - focus_start);
        ratio2_fwhm = zeros(41,focus_end - focus_start);
        ratio_fwhm = zeros(focus_end - focus_start,1);
        diameter_fwhm = zeros(focus_end - focus_start,1);
    end

    l=0;
    for k=focus_start:focus_end
        
        layer = processedImages(:,:,k);
        layer(layer<threshold) = 0;
        cog_layer = layer;
        cog_layer(cog_layer<cog_threshold) = 0;

        binaryImage = true(size(cog_layer));
        measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
        cog = measurements(1).WeightedCentroid;
        max_y = round(cog(1));
        max_x = round(cog(2));

        l=l+1;
        for i=max(max_x)
            profil1 = double(layer(i,:));
            rad1=1;
            max_profil1 = find(profil1 == max(profil1));
            [value, ind] = find(profil1>0.5*max(profil1));
            diameter1_fwhm(l) = length(ind);
            ratio1_fwhm(l) = max(profil1)/diameter1_fwhm(l);

        end
        j=0;
        for i=max(max_y)
            profil2 = double(layer(:,i));
            max_profil2 = find(profil2 == max(profil2));
            [ind, value] = find(profil2>0.5*max(profil2));
            diameter2_fwhm(l) = length(ind);
            ratio2_fwhm(l) = max(profil2)/diameter2_fwhm(l);

        end
        diameter_fwhm(l) = (diameter1_fwhm(l) + diameter2_fwhm(l))/2;
        ratio_fwhm(l) = double(max(max(layer))) / diameter_fwhm(l);
    end
    [v, focusz] = min(ratio_fwhm);
    val_fwhm = diameter_fwhm(focusz);
    fd_fwhm = round(val_fwhm*pixel_size,2);
    wd_fwhm = round((y_offset-focus_start-focusz+nozzle_offset)*pixel_size,2);

    disp(['Focus diameter FWHM: ',num2str(fd_fwhm), ' mm; Working distance FWHM: ', num2str(wd_fwhm), ' mm' ])

    focus_position = y_offset - (y_offset - focus_start - focusz);
    focus_layer = processedImages(:,:,focus_position);

end

cog_layer = focus_layer;
cog_layer(cog_layer<cog_threshold) = 0;
binaryImage = true(size(cog_layer));
measurements = regionprops(binaryImage, cog_layer, 'WeightedCentroid');
cog = measurements(1).WeightedCentroid;
max_y = round(cog(1));
max_x = round(cog(2));
figure
imshow(focus_layer, Colormap=jet)
viscircles([max_y max_x], round(val_fwhm/2));

%% x-y-focus position detection

binaryImage = true(size(focus_layer));
measurements = regionprops(binaryImage, focus_layer>threshold, 'WeightedCentroid');
cog = measurements(1).WeightedCentroid;
center_y = round(cog(1));
center_x = round(cog(2));

%% injector intensities for multijet nozzles

nozzle_position = y_offset + nozzle_offset;
focus_layer = processedImages(:,:,focus_position);

distance_25_percent = nozzle_position - ceil(0.25 * (nozzle_position-focus_position));
distance_25_percent_layer = processedImages(:,:,distance_25_percent);

if multijet_nozzle == 1
    distance_25_percent_layer_bw = distance_25_percent_layer;
    distance_25_percent_layer_bw(distance_25_percent_layer_bw<threshold_multijet)=0;
    distance_25_percent_layer_bw(distance_25_percent_layer_bw>threshold_multijet-1)=1;

    cc = bwconncomp(distance_25_percent_layer_bw, 8);
    L = labelmatrix(cc);
    rgb = label2rgb(L);

    numArea = unique(L);
    elements_of_area = zeros(length(numArea),2);
    for j=1:length(numArea)
        elements_of_area(j,1) = sum(sum(L==j));
        elements_of_area(j,2) = j;
    end
    sortiert = flip(sortrows(elements_of_area,1));
  
    intensity_injectors = zeros(number_injectors,1);
    text_injector = zeros(number_injectors, 3);
    texts = strings(number_injectors, 1);

    for i=1:number_injectors
        intensity_injectors(i) = sum(sum(distance_25_percent_layer(L==sortiert(i,2))));
    end
    total_intensity_injectors = sum(intensity_injectors);
    percent_injector = round(intensity_injectors/total_intensity_injectors,3);

    for i=1:number_injectors
        text_injector(i,1) = 255;
        [start_text1, start_text2, v] = find(L==sortiert(i,2));
        text_injector(i,3) = max(start_text1);
        text_injector(i,2) = floor((max(start_text2) + min(start_text2))/2);
        texts(i) = num2str(percent_injector(i));
    end
    layer_text = distance_25_percent_layer;
    layer_text = insertText(layer_text, text_injector(:, 2:3), texts, 'FontSize', 20, 'AnchorPoint', 'CenterTop', 'BoxColor', 'w');

    figure
    imshow(layer_text)

%     exportgraphics(gcf, 'Injector intensities.jpg')
end

%% figure for overview
x_coord1 = linspace(-x_offset*pixel_size, x_offset*pixel_size, size(processedImages,2));
y_coord1 = linspace(-size(processedImages,1)/2*pixel_size, size(processedImages,1)/2*pixel_size, size(processedImages,1));
x_coord2 = linspace(-size(processedImages,1)/2*pixel_size, size(processedImages,1)/2*pixel_size, size(processedImages,1));
y_coord2 = linspace(-nozzle_offset*pixel_size, -y_offset*pixel_size, y_offset);
x_coord3 = linspace(-x_offset*pixel_size, x_offset*pixel_size, size(processedImages,2));
y_coord3 = linspace(-nozzle_offset*pixel_size, -y_offset*pixel_size, size(processedImages,3));

[X1, Y1] = meshgrid(x_coord1, y_coord1);
[X2, Y2] = meshgrid(x_coord2, y_coord2);
[X3, Y3] = meshgrid(x_coord3, y_coord3);

figure('Position', [10 10 800 1000])
tiledlayout(2,2)
nexttile
surf(X1, Y1, double(flip(processedImages(:,:,focus_position))))
shading interp
set(gca, colormap=jet)
view(2)
axis equal
title('Focus Layer', ['z = ', num2str( round(-(y_offset-focus_position+nozzle_offset) * pixel_size, 2 )), ' mm'], FontSize=12)
xlim([min(x_coord1) max(x_coord1)])
ylim([min(y_coord1) max(y_coord1)])
xlabel('x / mm', 'Interpreter','latex')
ylabel('y / mm', 'Interpreter','latex')
clim([0 255])

nexttile
surf(X1, Y1, double(flip(distance_25_percent_layer)))
shading interp
set(gca, colormap=jet)
view(2)
axis equal
title('Layer 25 % Focus Position', ['z = ', num2str( round(-(y_offset-focus_position+nozzle_offset) * pixel_size/2, 2)), ' mm'], FontSize=12)
xlim([min(x_coord1) max(x_coord1)])
ylim([min(y_coord1) max(y_coord1)])
xlabel('x / mm', 'Interpreter','latex')
ylabel('y / mm', 'Interpreter','latex')
clim([0 255])

nexttile
surf(X3, Y3, flip(double(floor(xz_cross(:,:,center_x)))'));
shading interp
set(gca, colormap=jet)
view(2)
axis equal
title('x-z-Cross Section', ['y = ', num2str(round(y_coord1(size(y_coord1,2)-center_x),2)), ' mm'], FontSize=12)
xlim([min(x_coord3) max(x_coord3)])
ylim([min(y_coord3) max(y_coord3)])
xlabel('x / mm', 'Interpreter','latex')
ylabel('z / mm', 'Interpreter','latex')
clim([0 255])

nexttile
surf(X2, Y2, flip(double(floor(yz_cross(:,:,center_y)))'));
shading interp
set(gca, colormap=jet)
view(2)
axis equal
title('y-z-Cross Section', ['x = ', num2str(round(x_coord1(center_y),2)), ' mm'], FontSize=12)
xlim([min(x_coord2) max(x_coord2)])
ylim([min(y_coord2) max(y_coord2)])
xlabel('x / mm', 'Interpreter','latex')
ylabel('z / mm', 'Interpreter','latex')
clim([0 255])

cb = colorbar;
cb.Layout.Tile = 'south';

%% export graphics
% exportgraphics(gcf, 'Measurement overview.jpg')
