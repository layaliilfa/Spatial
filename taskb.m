data = readtable('SFbay_velocities.csv');
head(data)     % lihat 8 baris pertama

%% === STEP 1: BACA DAN SIAPKAN DATA ===
T = readtable('SFbay_velocities.csv');

Ve = T.Velocity_E_;
Vn = T.Velocity_N_;

figure;
scatter(Ve, Vn, 30, 'filled');
xlabel('Velocity East (mm/yr)');
ylabel('Velocity North (mm/yr)');
title('Velocity Space (Simpson et al., 2012)');
grid on; axis equal;


%% === TASK B: Cluster Analysis (Pure MATLAB Core) ===
clc; clear; close all;

% === 1. Load data ===
T = readtable('SFbay_velocities.csv');
Ve = T.Velocity_E_;
Vn = T.Velocity_N_;
X = [Ve, Vn];
X = double(X);

% === 2. Manual K-means Clustering (2–9 clusters) ===
maxK = 9;
colors = lines(maxK); % palet warna dasar
figure;
tiledlayout(2,4,'Padding','compact','TileSpacing','compact');

for k = 2:maxK
    % --- inisialisasi centroid acak
    rng(1); % biar hasil stabil
    C = X(randperm(size(X,1), k), :);
    
    for iter = 1:100
        % Hitung jarak tiap titik ke tiap centroid
        D = zeros(size(X,1), k);
        for j = 1:k
            D(:,j) = sqrt(sum((X - C(j,:)).^2, 2));
        end
        
        % Assign cluster terdekat
        [~, idx] = min(D, [], 2);
        
        % Update centroid
        newC = zeros(size(C));
        for j = 1:k
            newC(j,:) = mean(X(idx == j, :), 1);
        end
        
        % Stop jika konvergen
        if all(abs(newC - C) < 1e-6, 'all')
            break
        end
        C = newC;
    end
    
    % --- plot hasil cluster manual
    nexttile;
    hold on;
    for j = 1:k
        scatter(X(idx==j,1), X(idx==j,2), 20, colors(j,:), 'filled');
    end
    scatter(C(:,1), C(:,2), 50, 'k', 'p', 'filled'); % centroid hitam
    hold off;
    
    title(['k = ', num2str(k)]);
    xlabel('Velocity East (mm/yr)');
    ylabel('Velocity North (mm/yr)');
    axis equal; grid on;
end

sgtitle('Manual Cluster Analysis (Velocity Space, 2–9)');

%% === SIMPSON-STYLE CLUSTER FIGURE (FINAL MATCHED AXES) ===
figure('Position',[200 200 1200 500])

% === LEFT PANEL: Velocity Space ===
subplot(1,2,1)
hold on
for i = 1:k
    scatter(Ve(idx==i), Vn(idx==i), 30, cols(i,:), 'filled')
end
scatter(C(:,1), C(:,2), 60, 'k', 'p', 'filled')
xlabel('Velocity East (mm/yr)')
ylabel('Velocity North (mm/yr)')
title('GPS horizontal velocities (k = 4)')
axis equal; grid on

% >>> set axis limits to match Simpson (2012)
xlim([-30 5])
ylim([0 40])

% Tambahkan label blok
for i = 1:k
    text(mean(Ve(idx==i)), mean(Vn(idx==i)), ...
        blockNames{i}, 'FontWeight','bold','FontSize',9,...
        'HorizontalAlignment','center');
end
hold off

% === RIGHT PANEL: Spatial Map ===
subplot(1,2,2)
hold on
for i = 1:k
    scatter(Lon(idx==i), Lat(idx==i), 35, cols(i,:), 'filled')
end
xlabel('Longitude'); ylabel('Latitude');
title('Map of Locations (k = 4)')
grid on; axis equal

% >>> set map limits also consistent with paper
xlim([-123.5 -120])
ylim([36 40.5])
set(gca,'YDir','normal')

% Label blok di lokasi rata-rata
for i = 1:k
    text(mean(Lon(idx==i)), mean(Lat(idx==i)), ...
        extractBefore(blockNames{i},'('), ...
        'FontWeight','bold','FontSize',9,...
        'HorizontalAlignment','center');
end
hold off

sgtitle('Cluster Analysis for GNSS (Matched to Simpson et al., 2012)')

%% === TASK B: Velocity-Space Clustering Exploration (k = 2–10) ===
clc; clear; close all;

% === 1. Load GNSS data ===
T = readtable('SFbay_velocities.csv');
Ve = T.Velocity_E_;
Vn = T.Velocity_N_;
X = double([Ve, Vn]);

% === 2. Define plotting style ===
maxK = 10;
cols = lines(maxK);   % color palette

% === 3. Create figure layout ===
figure('Position',[100 100 1200 650],'Color','w');
tiledlayout(3,3,'Padding','compact','TileSpacing','compact');

% === 4. Manual K-means loop ===
for k = 2:maxK
    % Initialize centroids randomly
    rng(1);
    C = X(randperm(size(X,1), k), :);

    % Simple manual iteration (no toolbox)
    for iter = 1:100
        % Distance to centroids
        D = zeros(size(X,1), k);
        for j = 1:k
            D(:,j) = sqrt(sum((X - C(j,:)).^2, 2));
        end
        % Assign cluster
        [~, idx] = min(D, [], 2);
        % Update centroids
        newC = zeros(size(C));
        for j = 1:k
            newC(j,:) = mean(X(idx==j,:), 1);
        end
        if all(abs(newC - C) < 1e-6, 'all')
            break
        end
        C = newC;
    end

    % === Plot each k result ===
    nexttile;
    hold on
    for j = 1:k
        scatter(X(idx==j,1), X(idx==j,2), 25, cols(j,:), 'filled', ...
            'MarkerFaceAlpha',0.8, 'MarkerEdgeColor','none');
    end
    scatter(C(:,1), C(:,2), 50, 'k', 'p', 'filled'); % centroids
    hold off

    % Axis formatting
    title(['k = ', num2str(k)], 'FontWeight','bold', 'FontSize',11);
    xlabel('Velocity East (mm/yr)', 'FontSize',9);
    ylabel('Velocity North (mm/yr)', 'FontSize',9);
    xlim([-30 5]);
    ylim([0 40]);
    axis equal
    grid on
    set(gca,'FontSize',9,'LineWidth',0.8,'GridAlpha',0.25);
end

% === 5. Figure title ===
sgtitle({'Velocity-Space K-Means Clustering (k = 2–10)', ...
    'Scaled to Simpson et al. (2012) — Clean Version'}, ...
    'FontWeight','bold','FontSize',13);
%% === TASK B: Best-k via WCSS (Elbow) — no toolbox ===
clc; clear; close all;

% 1) Load data (Ve,Vn) — sama seperti sebelumnya
T  = readtable('SFbay_velocities.csv');
Ve = T.Velocity_E_;
Vn = T.Velocity_N_;
X  = double([Ve, Vn]);

% 2) Hitung WCSS untuk k = 2..10 dengan k-means manual
klist = 2:10;
WCSS  = zeros(size(klist));
idxAll = cell(size(klist));  % kalau mau pakai lagi

for ii = 1:numel(klist)
    k = klist(ii);
    [idx, C, wcss] = manual_kmeans_core(X, k, 100, 1); % maxIter=100, seed=1
    WCSS(ii)  = wcss;
    idxAll{ii}= idx;
end

% 3) Deteksi "elbow" otomatis (max distance to line k=min..max)
x1 = klist(1);   y1 = WCSS(1);
x2 = klist(end); y2 = WCSS(end);
% jarak tegak lurus tiap titik ke garis (x1,y1)->(x2,y2)
num   = abs((y2-y1).*klist - (x2-x1).*WCSS - x1*y2 + x2*y1);
den   = sqrt((y2-y1)^2 + (x2-x1)^2);
dist  = num./den;
[~,iElbow] = max(dist);
bestK = klist(iElbow);

% 4) Plot Elbow
figure('Color','w','Position',[200 200 700 420]);
plot(klist, WCSS, '-o', 'LineWidth',1.8, 'MarkerFaceColor',[0.2 0.5 1], 'MarkerSize',6);
hold on
plot(bestK, WCSS(iElbow), 'p', 'MarkerFaceColor','r', 'MarkerEdgeColor','k', 'MarkerSize',12);
yline(WCSS(iElbow), ':', 'Color',[0.4 0.4 0.4]);
xline(bestK, ':', 'Color',[0.4 0.4 0.4]);
title('Elbow Method using WCSS (manual k-means)','FontWeight','bold');
xlabel('Number of clusters, k'); ylabel('WCSS (sum of squared distances)');
grid on; box on
text(bestK+0.1, WCSS(iElbow), sprintf('  \\leftarrow best k = %d', bestK), ...
    'FontWeight','bold','Color','r');

% 5) (Opsional) Tabel ringkas pengurangan WCSS antar k
dWCSS = [NaN, diff(WCSS)];
relDrop = [NaN, -diff(WCSS)./WCSS(1:end-1)]; % proporsi penurunan
Tsummary = table(klist(:), WCSS(:), dWCSS(:), relDrop(:), ...
    'VariableNames', {'k','WCSS','Delta_WCSS','RelDrop'});
disp(Tsummary);

%% ====== local function: simple manual kmeans ======
function [idx, C, wcss] = manual_kmeans_core(X, k, maxIter, seed)
    if nargin < 3, maxIter = 100; end
    if nargin < 4, seed = 1; end
    rng(seed);

    % init centroid acak dari titik data
    C = X(randperm(size(X,1), k), :);

    for it = 1:maxIter
        % jarak Euclidean ke tiap centroid
        D = zeros(size(X,1), k);
        for j = 1:k
            diffj = X - C(j,:);
            D(:,j) = sqrt(sum(diffj.^2, 2));
        end
        % assign cluster terdekat
        [~, idx] = min(D, [], 2);

        % update centroid
        newC = zeros(k, size(X,2));
        for j = 1:k
            pts = X(idx==j, :);
            if isempty(pts)
                % jika ada cluster kosong, re-init centroid ke titik acak
                newC(j,:) = X(randi(size(X,1)), :);
            else
                newC(j,:) = mean(pts, 1);
            end
        end

        if all(abs(newC - C) < 1e-6, 'all'), break; end
        C = newC;
    end

    % hitung WCSS
    wcss = 0;
    for j = 1:k
        pts = X(idx==j,:);
        if ~isempty(pts)
            dif = pts - C(j,:);
            wcss = wcss + sum(sum(dif.^2));
        end
    end
end


%% === FINAL CLUSTER VISUALIZATION (k = 4, NO TOOLBOX, CLEAN VERSION) ===
clc; clear; close all;

% === 1. Load data ===
T = readtable('SFbay_velocities.csv');
Ve = T.Velocity_E_;
Vn = T.Velocity_N_;
Lon = T.Longitude;
Lat = T.Latitude;
X = [Ve, Vn];

% === 2. Manual k-means ===
k = 4;                    % jumlah klaster optimal
rng(1);                   % seed random biar konsisten
C = X(randperm(size(X,1), k), :);  % inisialisasi centroid acak

for iter = 1:100
    % Hitung jarak manual (tanpa pdist2)
    D = zeros(size(X,1), k);
    for j = 1:k
        D(:,j) = sqrt( (X(:,1)-C(j,1)).^2 + (X(:,2)-C(j,2)).^2 );
    end

    % Assign cluster terdekat
    [~, idx] = min(D, [], 2);

    % Update centroid
    newC = zeros(size(C));
    for j = 1:k
        newC(j,:) = mean(X(idx==j,:), 1);
    end

    % Berhenti jika sudah stabil
    if all(abs(newC - C) < 1e-6, 'all')
        break
    end
    C = newC;
end

% === 3. Definisikan warna & marker ===
cols = lines(k);
markers = {'^','s','o','d'};  % segitiga, persegi, bulat, belah ketupat
blockNames = {'Pacific Block','Bay Block','East Bay Block','Sierra Nevada–Great Valley Block'};

%% === 4. Velocity-Space Plot (skala Simpson et al. 2012) ===
figure('Position',[200 200 700 550],'Color','w');
hold on;
for i = 1:k
    % Pastikan idx ada titik untuk klaster ini
    if any(idx == i)
        scatter(Ve(idx==i), Vn(idx==i), 40, cols(i,:), markers{i}, ...
            'filled','MarkerEdgeColor','k');
    end
end
% Plot centroid
scatter(C(:,1), C(:,2), 60, 'k', 'p', 'filled');

xlabel('Velocity East (mm/yr)');
ylabel('Velocity North (mm/yr)');
title('Velocity-Space Clusters (k = 4) — Scaled to Simpson et al., 2012','FontWeight','bold');
xlim([-5 30]); ylim([0 40]); axis equal;
grid on; box on;

for i = 1:k
    text(mean(Ve(idx==i)), mean(Vn(idx==i)), blockNames{i}, ...
        'FontWeight','bold','FontSize',9,'HorizontalAlignment','center');
end
legend(blockNames,'Location','southeastoutside','Box','off');
hold off;

%% === 5. Spatial Map (Longitude–Latitude) ===
figure('Position',[200 200 800 600],'Color','w');
hold on;

% Garis San Andreas (dummy)
sanAndreasLat = [36.5 37 37.5 38 38.5 39 39.5 40];
sanAndreasLon = [-122.5 -122.4 -122.3 -122.4 -122.35 -122.4 -122.5 -122.45];
plot(sanAndreasLon, sanAndreasLat, 'r-', 'LineWidth',1.5);

% Titik GNSS tiap klaster
for i = 1:k
    if any(idx == i)
        scatter(Lon(idx==i), Lat(idx==i), 40, cols(i,:), markers{i}, ...
            'filled','MarkerEdgeColor','k');
    end
end

% Label blok di tengah-tengah titik
for i = 1:k
    text(mean(Lon(idx==i)), mean(Lat(idx==i)), blockNames{i}, ...
        'FontWeight','bold','FontSize',9,'HorizontalAlignment','center');
end

xlabel('Longitude'); ylabel('Latitude');
xlim([-123.5 -120]); ylim([36 40.5]); axis equal;
title('GNSS Cluster Map (k = 4) — San Andreas Fault Region','FontWeight','bold');
legend(blockNames,'Location','southeastoutside','Box','off');
grid on;

%% === MAP OF LOCATIONS (N clusters = 4) — Refined realistic geometry ===
clc; clear; close all;

% 1. Load data
T  = readtable('SFbay_velocities.csv');
Lon = T.Longitude;
Lat = T.Latitude;
Ve  = T.Velocity_E_;
Vn  = T.Velocity_N_;
X   = [Ve Vn];

% 2. Manual K-means
k = 4; rng(1);
C = X(randperm(size(X,1),k),:);
for it = 1:100
    D = zeros(size(X,1),k);
    for j=1:k
        D(:,j) = sqrt((X(:,1)-C(j,1)).^2 + (X(:,2)-C(j,2)).^2);
    end
    [~,idx] = min(D,[],2);
    newC = zeros(size(C));
    for j=1:k, newC(j,:) = mean(X(idx==j,:),1); end
    if all(abs(newC-C)<1e-6,'all'), break; end
    C = newC;
end

% 3. Warna & marker
cols = [0.8 0.1 0.1; 0.95 0.8 0.1; 0.1 0.8 0.2; 0.55 0.35 0.8]; 
markers = {'s','^','d','o'};
blockNames = {'PB','BB','EBB','SNGVB'};

% 4. Fault lines (koordinat lebih akurat & halus)
sanAndreasLon = [-123.45 -123.35 -123.25 -123.10 -122.95 -122.8 -122.65 -122.55 -122.45 -122.35 -122.25 -122.1];
sanAndreasLat = [36.3 36.6 36.8 37.1 37.4 37.7 38.0 38.4 38.8 39.1 39.3 39.5];

haywardLon = [-122.05 -121.98 -121.92 -121.85 -121.78 -121.72];
haywardLat = [37.2 37.5 37.9 38.3 38.7 39.1];

calaverasLon = [-121.8 -121.72 -121.65 -121.58 -121.50];
calaverasLat = [36.7 37.1 37.6 38.2 38.7];

% 5. Outline pantai sederhana
coastLon = [-123.5 -123.3 -123.1 -122.9 -122.7 -122.6 -122.5 -122.4 -122.3 -122.2 -122.2];
coastLat = [36.3 36.6 36.9 37.2 37.4 37.6 37.8 38.1 38.4 38.8 39.2];

% 6. Plot
figure('Position',[100 100 750 900],'Color','w');
hold on;

% Coastline (abu)
plot(coastLon, coastLat, 'Color',[0.4 0.4 0.4],'LineWidth',0.8);

% Fault lines hijau lembut
plot(sanAndreasLon, sanAndreasLat, 'Color',[0 0.6 0],'LineWidth',1.3);
plot(haywardLon, haywardLat, 'Color',[0 0.6 0],'LineWidth',1.3);
plot(calaverasLon, calaverasLat, 'Color',[0 0.6 0],'LineWidth',1.3);

% Plot titik cluster
order = [4 3 2 1];
for i = order
    scatter(Lon(idx==i), Lat(idx==i), 45, cols(i,:), markers{i}, ...
        'filled','MarkerEdgeColor','k','LineWidth',0.6);
end

% Label blok (digeser biar rapi)
labelShift = [-0.07 0.02; -0.07 0.07; -0.05 0.05; 0 0];
for i = 1:k
    text(mean(Lon(idx==i))+labelShift(i,1), mean(Lat(idx==i))+labelShift(i,2), ...
        blockNames{i}, 'FontWeight','bold','FontSize',10, ...
        'HorizontalAlignment','center','Color','k');
end

xlabel('Longitude','FontWeight','bold');
ylabel('Latitude','FontWeight','bold');
title('Map of Locations (N clusters = 4)','FontWeight','bold','FontSize',12);

xlim([-123.5 -120]);
ylim([36 40.5]);
set(gca, 'XDir', 'reverse');   % <=== ini yang penting biar barat di kanan
axis equal; 
pbaspect([1 1.35 1]);
grid on; box on;

set(gca,'FontSize',10,'LineWidth',1,'XMinorTick','on','YMinorTick','on');
