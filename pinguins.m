data = readtable('penguins_.csv');
disp(data.Properties.VariableNames);

% Cek jumlah missing values di tiap kolom
sum(ismissing(data))

% Buat boxplot
figure;
boxplot([data.flipper_length_mm, data.body_mass_g], ...
        'Labels', {'Flipper length','Body mass'});
title('Boxplot Flipper Length vs Body Mass');
ylabel('Value');
grid on;

% --- Ambil data ---
flipper = data.flipper_length_mm;
body    = data.body_mass_g;

% --- Hapus missing values dulu ---
flipper = flipper(~ismissing(flipper));
body    = body(~ismissing(body));

% --- Hitung Q1, Q3, dan IQR ---
Q1_fl = quantile(flipper, 0.25);
Q3_fl = quantile(flipper, 0.75);
IQR_fl = Q3_fl - Q1_fl;

Q1_bd = quantile(body, 0.25);
Q3_bd = quantile(body, 0.75);
IQR_bd = Q3_bd - Q1_bd;

% --- Tentukan batas bawah & atas ---
lower_fl = Q1_fl - 1.5 * IQR_fl;
upper_fl = Q3_fl + 1.5 * IQR_fl;

lower_bd = Q1_bd - 1.5 * IQR_bd;
upper_bd = Q3_bd + 1.5 * IQR_bd;

% --- Hitung jumlah outlier ---
outlier_fl = sum(flipper < lower_fl | flipper > upper_fl);
outlier_bd = sum(body < lower_bd | body > upper_bd);

fprintf('Jumlah outlier Flipper length: %d\n', outlier_fl);
fprintf('Jumlah outlier Body mass: %d\n', outlier_bd);

figure;
boxplot([flipper, body], 'Labels', {'Flipper length','Body mass'}, ...
    'Whisker', 1.5, ...             % gunakan standar 1.5 IQR
    'Symbol', 'r+', ...             % warna merah untuk outlier
    'Colors', ['b','k']);           % warna kotak biru & hitam
title('Boxplot of Flipper Length and Body Mass (with Outliers)');
ylabel('Value');
grid on;

% Tambah garis median (opsional)
hold on;
medians = [median(flipper), median(body)];
plot(1:2, medians, 'go', 'MarkerSize', 8, 'MarkerFaceColor', 'g');

% === STEP 2: REMOVE OUTLIERS AND NON-IDENTIFIED DATA ===

% Hapus data tidak teridentifikasi (missing values)
data_clean = rmmissing(data);

% Ambil dua variabel utama
flipper = data_clean.flipper_length_mm;
body    = data_clean.body_mass_g;

% Hitung Q1, Q3, dan IQR
Q1_fl = quantile(flipper, 0.25);
Q3_fl = quantile(flipper, 0.75);
IQR_fl = Q3_fl - Q1_fl;

Q1_bd = quantile(body, 0.25);
Q3_bd = quantile(body, 0.75);
IQR_bd = Q3_bd - Q1_bd;

% Tentukan batas bawah dan atas
lower_fl = Q1_fl - 1.5 * IQR_fl;
upper_fl = Q3_fl + 1.5 * IQR_fl;
lower_bd = Q1_bd - 1.5 * IQR_bd;
upper_bd = Q3_bd + 1.5 * IQR_bd;

% Hapus data di luar rentang IQR
validRows = (flipper >= lower_fl & flipper <= upper_fl) & ...
            (body >= lower_bd & body <= upper_bd);
data_final = data_clean(validRows, :);

fprintf('Jumlah data awal: %d\n', height(data));
fprintf('Jumlah data setelah dibersihkan: %d\n', height(data_final));

% === BOX PLOT BEFORE & AFTER CLEANING ===

% Data awal (belum dibersihkan)
flipper_raw = data.flipper_length_mm;
body_raw    = data.body_mass_g;

% Data bersih (setelah hapus missing & outlier)
flipper_clean = data_final.flipper_length_mm;
body_clean    = data_final.body_mass_g;

% Gabungkan untuk plotting
figure;

subplot(1,2,1)
boxplot([flipper_raw, body_raw], 'Labels', {'Flipper length','Body mass'}, ...
    'Symbol','r+');
title('Before Cleaning');
ylabel('Value');
grid on;

subplot(1,2,2)
boxplot([flipper_clean, body_clean], 'Labels', {'Flipper length','Body mass'}, ...
    'Symbol','r+');
title('After Cleaning');
grid on;

% === STEP 3: K-MEANS CLUSTERING ANALYSIS ===

% Ambil dua variabel utama dari data bersih
X = [data_final.flipper_length_mm, data_final.body_mass_g];

% Inisialisasi
maxK = 10;              % rentang nilai k
WCCS = zeros(maxK,1);   % untuk menyimpan total sum of squares tiap k

% Loop untuk menghitung WCCS setiap k
for k = 1:maxK
    [~, ~, sumd] = kmeans(X, k, 'Replicates',5, 'Display','off');
    WCCS(k) = sum(sumd);  % jumlahkan seluruh cluster sum of squares
end

% Plot Elbow Method
figure;
plot(1:maxK, WCCS, '-o', 'LineWidth',1.5);
xlabel('Number of Clusters (k)');
ylabel('Within-Cluster Sum of Squares (WCCS)');
title('Elbow Method to Determine Optimum k');
grid on;

% === STEP 4: CLUSTERING WITH OPTIMUM k ===
k_opt = 3; % hasil dari Elbow Method

% Jalankan K-Means clustering
[idx, C] = kmeans(X, k_opt, 'Replicates', 5, 'Display', 'final');

% Scatter plot hasil klasterisasi
figure;
gscatter(X(:,1), X(:,2), idx, 'rgb', 'o', 8);
hold on;
plot(C(:,1), C(:,2), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
xlabel('Flipper Length (mm)');
ylabel('Body Mass (g)');
title(sprintf('K-Means Clustering Result (k = %d)', k_opt));
legend('Cluster 1', 'Cluster 2', 'Cluster 3', 'Centroids');
grid on;
