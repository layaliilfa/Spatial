%% ============================================================
%  TASK A – IMAGE CLUSTERING USING K-MEANS
%  Nomor 2: Data Preparation (konversi RGB & cleaning)
%  Nomor 3: Elbow Method untuk menentukan k optimum
% ============================================================

clc; clear; close all;

%% === STEP 1: BACA DAN TAMPILKAN CITRA ===
% Hindari konflik nama variabel 'img' dengan file .m
rgbImg = imread('img.jpg');
figure;
imshow(rgbImg);
title('Original Image');

% Konversi ke double precision (range 0–1)
img_double = im2double(rgbImg);

% Dapatkan ukuran citra
[nRows, nCols, nBands] = size(img_double);
fprintf('Ukuran citra: %d x %d pixels\n', nRows, nCols);

%% === STEP 2: REKONSTRUKSI DATA RGB (n x 3) ===
% Setiap baris merepresentasikan 1 pixel, dan setiap kolom 3 kanal warna (R,G,B)
pixelData = reshape(img_double, nRows*nCols, nBands);

% Jika citra besar, sampling sebagian pixel biar lebih cepat
sampleSize = 50000; % ambil 50 ribu pixel acak
sample_idx = randperm(size(pixelData,1), sampleSize);
pixelSample = pixelData(sample_idx, :);

%% === STEP 3: TENTUKAN NILAI k OPTIMUM (ELBOW METHOD) ===
maxK = 10;
WCCS = zeros(maxK,1);

for k = 1:maxK
    [idx, C, sumd] = kmeans(pixelSample, k, ...
        'Replicates', 5, 'MaxIter', 300, 'Display', 'off');
    WCCS(k) = sum(sumd);
end

% Plot Elbow Curve
figure;
plot(1:maxK, WCCS, '-o', 'LineWidth', 1.5);
xlabel('Number of Clusters (k)');
ylabel('Within-Cluster Sum of Squares (WCCS)');
title('Elbow Method for Image Clustering');
grid on;

%% === STEP 4: CLUSTERING PENUH DENGAN NILAI k OPTIMUM ===
% Misal dari grafik terlihat elbow di k = 4 (ubah sesuai hasil kamu)
k_opt = 4;

[idx_full, C_full] = kmeans(pixelData, k_opt, ...
    'Replicates', 5, 'MaxIter', 300, 'Display', 'off');

% Rekonstruksi hasil klaster jadi citra baru
segmented = reshape(C_full(idx_full, :), nRows, nCols, nBands);

figure;
imshow(segmented);
title(sprintf('Segmented Image (k = %d)', k_opt));

%% === STEP 5: INFORMASI TAMBAHAN ===
fprintf('Jumlah pixel: %d\n', nRows*nCols);
fprintf('Jumlah klaster (k): %d\n', k_opt);
disp('Centroid warna tiap klaster (RGB):');
disp(C_full);

%% === VISUALISASI CENTROID WARNA ===
figure;
colorBar = uint8(reshape(C_full * 255, [1, size(C_full,1), 3]));  % ubah ke format RGB 0–255
imshow(colorBar);
title('Centroid Colors (RGB of Each Cluster)');
xlabel('Cluster Index');

%% === MEMBUAT PETA TUTUPAN LAHAN (KATEGORIK) ===

% Buat label kelas manual sesuai interpretasi klaster
% (sesuaikan urutan klaster dengan warna centroid yang kamu lihat)
% Misal: 1=Hutan, 2=Permukiman, 3=Vegetasi sedang, 4=Air/Bayangan
labels = ["Hutan", "Permukiman", "Vegetasi sedang", "Air/Bayangan"];

% Warna representatif untuk tiap kelas (dalam RGB 0–1)
classColors = [
    0.0, 0.5, 0.0;   % hijau (hutan)
    0.8, 0.6, 0.4;   % coklat terang (permukiman)
    0.3, 0.8, 0.3;   % hijau muda (vegetasi sedang)
    0.1, 0.1, 0.3    % biru tua (air/bayangan)
];

% Buat image baru dengan warna berdasarkan label klaster
landCoverRGB = reshape(classColors(idx_full,:), nRows, nCols, 3);

figure;
imshow(landCoverRGB);
title('Peta Tutupan Lahan Berdasarkan K-Means');

% Tambah legenda warna di bawah
figure;
imshow(landCoverRGB);
title('Peta Tutupan Lahan');

% Buat bar warna (legend)
hold on;
for i = 1:length(labels)
    patch([0.02 0.07 0.07 0.02], [0.95-(i-1)*0.05 0.95-(i-1)*0.05 0.9-(i-1)*0.05 0.9-(i-1)*0.05], classColors(i,:), 'EdgeColor','none');
    text(0.08, 0.91-(i-1)*0.05, labels(i), 'Color','w', 'FontSize',10, 'FontWeight','bold');
end
hold off;

%% === STEP 6: PERBANDINGAN HASIL DENGAN k = 6 ===
k_compare = 6;

fprintf('\n=== Clustering tambahan untuk k = %d ===\n', k_compare);
[idx_k6, C_k6] = kmeans(pixelData, k_compare, ...
    'Replicates', 5, 'MaxIter', 300, 'Display', 'off');

segmented_k6 = reshape(C_k6(idx_k6, :), nRows, nCols, nBands);

% tampilkan hasil berdampingan
figure;
subplot(1,2,1);
imshow(segmented);
title(sprintf('Segmented Image (k = %d)', k_opt));

subplot(1,2,2);
imshow(segmented_k6);
title(sprintf('Segmented Image (k = %d)', k_compare));
sgtitle('Perbandingan Hasil Segmentasi (Elbow vs Visual)');

%% === CETAK CENTROID WARNA UNTUK k = 6 ===
fprintf('\nCentroid warna tiap klaster (RGB) untuk k = %d:\n', k_compare);
disp(C_k6);

%% === OPSIONAL: BUAT PETA TUTUPAN LAHAN DARI HASIL k=6 ===
% kamu boleh sesuaikan warna dan label sesuai interpretasi k=6
labels6 = ["Hutan lebat", "Vegetasi ringan", "Tanah terbuka", "Permukiman", "Air/Bayangan", "Area campuran"];

classColors6 = [
    0.0, 0.5, 0.0;   % hijau tua
    0.5, 0.8, 0.4;   % hijau muda
    0.9, 0.7, 0.4;   % coklat muda
    0.8, 0.4, 0.2;   % merah bata
    0.1, 0.1, 0.3;   % biru tua
    0.5, 0.5, 0.5    % abu campuran
];

landCoverRGB_k6 = reshape(classColors6(idx_k6,:), nRows, nCols, 3);

figure;
imshow(landCoverRGB_k6);
title('Peta Tutupan Lahan (k = 6)');

% buat legenda warna
hold on;
for i = 1:length(labels6)
    patch([0.02 0.07 0.07 0.02], [0.95-(i-1)*0.05 0.95-(i-1)*0.05 0.9-(i-1)*0.05 0.9-(i-1)*0.05], classColors6(i,:), 'EdgeColor','none');
    text(0.08, 0.91-(i-1)*0.05, labels6(i), 'Color','w', 'FontSize',10, 'FontWeight','bold');
end
hold off;
