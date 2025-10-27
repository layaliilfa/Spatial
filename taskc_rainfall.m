T = readtable('hujan.csv');
rain = T.mean_mm;
date = datetime(T.date);

meanRain = mean(rain, 'omitnan');
stdRain = std(rain, 'omitnan');
minRain = min(rain);
maxRain = max(rain);
medianRain = median(rain);

figure;
plot(date, rain, 'b-', 'LineWidth', 1.5);
xlabel('Tanggal'); ylabel('Curah Hujan (mm/hari)');
title('Variabilitas Harian Curah Hujan - Babakan Madang (Jul–Nov 2022)');
grid on;

% Histogram
figure;
histogram(rain, 30, 'FaceColor', [0.2 0.5 0.9]);
xlabel('Curah Hujan (mm/hari)');
ylabel('Frekuensi');
title('Distribusi Curah Hujan Harian');

fprintf('Rata-rata: %.2f mm/hari\n', meanRain);
fprintf('Simpangan baku: %.2f mm/hari\n', stdRain);
fprintf('Min–Max: %.2f – %.2f mm/hari\n', minRain, maxRain);
fprintf('Median: %.2f mm/hari\n', medianRain);


%% ============================================================
%  TASK C - Point 2
%  Empirical (Temporal) Covariance Function
%  Dataset: CHIRPS_BabakanMadang_Daily_Jul_Nov_2022.csv
%  ============================================================

clear; clc; close all;

%% 1) BACA DATA
% pastikan file CSV ada di folder kerja (Current Folder)
T = readtable('hujan.csv');

% lihat nama kolom kalau perlu
disp(T.Properties.VariableNames)

% ambil kolom tanggal & curah hujan
rain = T.mean_mm;                                   % kolom nilai curah hujan (mm/hari)
date = datetime(T.date,'InputFormat','yyyy-MM-dd'); % ubah jadi tipe tanggal

% tampilkan ringkasan cepat
fprintf('Jumlah data: %d hari, rentang %s s.d. %s\n', ...
        length(rain), string(min(date)), string(max(date)));

%% 2) CEK DAN BERSIHKAN NILAI HILANG
nanCount = sum(isnan(rain));
fprintf('Jumlah data hilang (NaN): %d\n', nanCount);
rain = fillmissing(rain,'linear');   % isi NaN secara linier kalau ada

%% 3) HITUNG EMPIRICAL (TEMPORAL) COVARIANCE
maxLag = 30;        % jendela lag (hari)
lags = 0:maxLag;
covVals = zeros(size(lags));

for k = 1:length(lags)
    lag = lags(k);
    if lag == 0
        covVals(k) = var(rain,'omitnan');  % kovarians nol-lag = varians
    else
        C = cov(rain(1:end-lag), rain(1+lag:end), 'omitrows');
        covVals(k) = C(1,2);               % ambil elemen kovarians antar-lag
    end
end

%% 4) NORMALISASI
covNorm = covVals / covVals(1);   % supaya Cov(0) = 1

%% 5) PLOT HASIL
figure;
stem(lags, covNorm, 'filled','LineWidth',1.3);
xlabel('Lag (hari)');
ylabel('Normalized Covariance');
title('Empirical Temporal Covariance Function - CHIRPS Babakan Madang (Jul–Nov 2022)');
grid on;

%% 6) INTERPRETASI OTOMATIS (SEDERHANA)
[~,decayIdx] = min(abs(covNorm - 0.1)); % lag saat covariance ≈ 0.1
fprintf('Korelasi menurun <0.1 pada lag ke-%d hari.\n',lags(decayIdx));

%% 7) (OPSIONAL) SIMPAN KE CSV
outTable = table(lags', covVals', covNorm', ...
                 'VariableNames',{'Lag_day','Covariance','NormalizedCov'});
writetable(outTable, 'Empirical_Covariance_CHIRPS_BabakanMadang.csv');
fprintf('File output disimpan sebagai Empirical_Covariance_CHIRPS_BabakanMadang.csv\n');


%% ============================================================
%  TASK C - Point 3
%  Fit Fungsi Kovarians Empiris dengan Model Analitik (Hirvonen)
%  ============================================================

clear; clc; close all;

%% 1) BACA DATA EMPIRIS
T = readtable('Empirical_Covariance_CHIRPS_BabakanMadang.csv');
lag = T.Lag_day;           % kolom lag (hari)
covEmp = T.NormalizedCov;  % kovarians empiris ternormalisasi

%% 2) MODEL ANALITIK: FUNGSI HIRVONEN / EKSPONENSIAL
% Bentuk model: C(h) = C0 * exp(-h / a)
hirvonenFun = @(p, h) p(1) * exp(-h ./ p(2));

% Parameter awal [C0, a]
initParams = [1, 5];

% Fungsi objektif (error kuadrat)
objFun = @(p) sum((hirvonenFun(p, lag) - covEmp).^2);

% Optimasi dengan fminsearch (tanpa toolbox tambahan)
fitParams = fminsearch(objFun, initParams);

C0_fit = fitParams(1);
a_fit  = fitParams(2);

fprintf('=== HASIL PARAMETER FIT ===\n');
fprintf('C0 = %.3f\n', C0_fit);
fprintf('a  = %.3f hari (range parameter)\n', a_fit);

%% 3) HITUNG NILAI FIT DAN RMSE
covModel = hirvonenFun(fitParams, lag);
RMSE = sqrt(mean((covModel - covEmp).^2,'omitnan'));
fprintf('RMSE = %.4f\n', RMSE);

%% 4) PLOT HASIL PERBANDINGAN
figure;
stem(lag, covEmp, 'b','filled','DisplayName','Empiris'); hold on;
plot(lag, covModel, 'r-','LineWidth',2,'DisplayName','Model Hirvonen');
xlabel('Lag (hari)');
ylabel('Normalized Covariance');
title('Fungsi Kovarians Empiris vs Analitik (Hirvonen) - CHIRPS Babakan Madang');
legend('Location','northeast');
grid on;

%% 5) TITIK DIMANA KONTRIBUSI SUDAH HILANG
idx_zero = find(covModel < 0.05, 1);
if ~isempty(idx_zero)
    fprintf('Kovarians < 0.05 pada lag ≈ %d hari (data tidak lagi saling berkontribusi)\n', lag(idx_zero));
else
    fprintf('Kovarians belum turun < 0.05 dalam jendela ini.\n');
end

%% 6) SIMPAN HASIL FIT
outTable = table(lag, covEmp, covModel, ...
    'VariableNames', {'Lag_day','EmpiricalCov','AnalyticalCov'});
writetable(outTable, 'Fitted_Hirvonen_Covariance.csv');
fprintf('File hasil fit disimpan sebagai Fitted_Hirvonen_Covariance.csv\n');

%% ============================================================
% TASK C - Poin 4 (versi realistis, hasil GEE 2 km)
% Surface curah hujan (grid 250 m) + Standard Error of Estimates (SEE)
% Input: CHIRPS_BabakanMadang_Mean_JulNov2022_Resampled2km.csv
% ============================================================

clear; clc; close all;

%% 1) BACA DATA TITIK HUJAN
T = readtable('CHIRPS_BabakanMadang_Mean_JulNov2022_Resampled2km.csv');
lon = T.lon;  lat = T.lat;  z = T.mean_mmday;
valid = ~isnan(lon) & ~isnan(lat) & ~isnan(z);
lon = lon(valid); lat = lat(valid); z = z(valid);
N = numel(z);
fprintf('Jumlah titik valid: %d\n', N);

%% 2) PROYEKSI KE METER (local tangent plane)
lon0 = mean(lon); lat0 = mean(lat);
R = 6371000;
x = deg2rad(lon - lon0) .* R .* cosd(lat0);
y = deg2rad(lat - lat0) .* R;

%% 3) GRID 250 m DI AOI
pad = 250;
xmin = min(x)-pad; xmax = max(x)+pad;
ymin = min(y)-pad; ymax = max(y)+pad;
dx = 250; dy = 250;
[xg, yg] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
[nr, nc] = size(xg);
fprintf('Ukuran grid 250 m: %d x %d sel\n', nr, nc);

%% 4) MODEL KOVARIANSI EKSPONENSIAL (Hirvonen-like, tanpa toolbox)
% Hitung jarak antar titik
D = zeros(N,N);
for i = 1:N
    for j = i:N
        d = sqrt((x(i)-x(j))^2 + (y(i)-y(j))^2);
        D(i,j) = d; D(j,i) = d;
    end
end

sigma2 = var(z,'omitnan');  % variansi global
a = 3000;                   % range spasial (3 km)
nugget = 0.05 * sigma2;     % noise kecil
fprintf('Parameter kovarian spasial: sigma2=%.3f, a=%.0f m, nugget=%.3f\n', sigma2, a, nugget);

K = sigma2 * exp(-D./a) + nugget * eye(N);
K = (K + K.')/2 + 1e-8 * eye(N);
L = chol(K,'lower');
alpha = L'\(L\z);

%% 5) INTERPOLASI GRID + SEE
Zg = zeros(nr,nc);
SEg = zeros(nr,nc);

for i = 1:nr
    for j = 1:nc
        % jarak dari titik ke grid node (i,j)
        d = sqrt((x - xg(i,j)).^2 + (y - yg(i,j)).^2);
        kstar = sigma2 * exp(-d./a);

        mu = kstar' * alpha;
        v  = L \ kstar;
        varp = sigma2 - sum(v.^2);
        if varp < 0, varp = 0; end

        Zg(i,j)  = mu;
        SEg(i,j) = sqrt(varp);
    end
end

%% 6) PLOT HASIL
figure;
subplot(1,2,1);
imagesc(xmin:xmax, ymin:ymax, Zg);
axis image xy; colorbar;
title('Surface Curah Hujan (mm/hari) – Grid 250 m');
xlabel('x (m)'); ylabel('y (m)');

subplot(1,2,2);
imagesc(xmin:xmax, ymin:ymax, SEg);
axis image xy; colorbar;
title('Standard Error of Estimates (mm/hari)');
xlabel('x (m)'); ylabel('y (m)');

%% 7) STATISTIK SCATTERED vs GRIDDED
stats = @(v)[mean(v,'omitnan') std(v,'omitnan') min(v) max(v)];
S_s = stats(z);   S_g = stats(Zg(:));
fprintf('\nStatistik Hujan (mm/hari)\n');
fprintf('Scattered (titik): mean=%.2f std=%.2f min=%.2f max=%.2f\n', S_s);
fprintf('Gridded (250 m):   mean=%.2f std=%.2f min=%.2f max=%.2f\n', S_g);

% Uji beda mean manual (tanpa toolbox)
m1 = mean(z,'omitnan');  s1 = std(z,'omitnan'); n1 = numel(z);
m2 = mean(Zg(:),'omitnan'); s2 = std(Zg(:),'omitnan'); n2 = numel(Zg);
tval = (m1 - m2) / sqrt((s1^2/n1) + (s2^2/n2));
pval = 2 * (1 - erf(abs(tval)/sqrt(2)));  % aproksimasi distribusi normal
pval = min(pval, 1);
fprintf('Uji beda mean (p≈%.4f) %s\n', pval, ternary(pval<0.05,'→ berbeda signifikan','→ tidak berbeda signifikan'));

%% 8) SIMPAN CSV OUTPUT
[Xout,Yout] = meshgrid(xmin:dx:xmax, ymin:dy:ymax);
out = table(Xout(:), Yout(:), Zg(:), SEg(:), 'VariableNames',{'x_m','y_m','rain_mmday','SE_mmday'});
writetable(out,'Surface_SEE_250m_BabakanMadang_Resampled2km.csv');
fprintf('File disimpan: Surface_SEE_250m_BabakanMadang_Resampled2km.csv\n');

%% === Helper fungsi ternary ===
function s = ternary(cond,a,b)
  if cond, s=a; else, s=b; end
end
