%% ============================
%  TASK C - POINT 1
%  Variability of Subsidence in Bandung
%  ============================

clear; clc; close all;

%% 1. Load data
load('data.mat');  % pastikan file data.mat ada di folder kerja

% Ambil field dari struktur
X = data.east;      % koordinat Easting (m)
Y = data.north;     % koordinat Northing (m)
Z = data.DV;        % nilai subsidence (mm atau cm)
Zs = data.DVs;      % standard deviation atau uncertainty

%% 2. Plot variabilitas subsidence
figure('Color','w');
scatter(X, Y, 40, Z, 'filled');
colorbar;
xlabel('Easting (m)', 'FontWeight', 'bold');
ylabel('Northing (m)', 'FontWeight', 'bold');
title('Subsidence Variability in Bandung', 'FontWeight', 'bold');
axis equal; grid on;

%% 3. Statistik deskriptif
meanZ   = mean(Z);
stdZ    = std(Z);
minZ    = min(Z);
maxZ    = max(Z);
medianZ = median(Z);

fprintf('\n=== Statistik Subsidence Bandung ===\n');
fprintf('Mean Subsidence      : %.3f\n', meanZ);
fprintf('Standard Deviation   : %.3f\n', stdZ);
fprintf('Minimum Subsidence   : %.3f\n', minZ);
fprintf('Maximum Subsidence   : %.3f\n', maxZ);
fprintf('Median Subsidence    : %.3f\n', medianZ);

%% 4. Titik dengan penurunan terbesar
[~, idx_max] = max(Z);
max_point = [X(idx_max), Y(idx_max), Z(idx_max)];
fprintf('\nTitik dengan subsidence terbesar:\n');
fprintf('E = %.2f m, N = %.2f m, Subsidence = %.3f\n', ...
        max_point(1), max_point(2), max_point(3));

% Tandai di plot
hold on;
plot(X(idx_max), Y(idx_max), 'kp', 'MarkerSize', 12, ...
     'MarkerFaceColor', 'r', 'DisplayName', 'Maximum Subsidence');
legend('Data Points', 'Maximum Subsidence', 'Location', 'best');

%% ==========================================================
%  TASK C - POINT 2 : EMPIRICAL COVARIANCE FUNCTION (ECF)
%  Replikasi tampilan seperti materi dosen ITB
%  Menampilkan:
%     (1) Scatter seluruh pasangan
%     (2) Garis mean tiap window (ECF halus)
%  ==========================================================

clear; clc; close all;
load('data.mat');  % berisi struct 'data' dengan field: east, north, DV

%% --- 1. Ambil data & normalisasi skala ---------------------------------
X = data.east;
Y = data.north;
Z = data.DV;

% Buang NaN
v = isfinite(X) & isfinite(Y) & isfinite(Z);
X = X(v); Y = Y(v); Z = Z(v);

% Jika standar deviasi kecil (asumsi satuan meter), ubah ke mm
if std(Z) < 0.1
    Z = Z * 1000;   % konversi meter → mm
end

% Zero-mean (tanpa detrending dulu)
Zc = Z - mean(Z);

fprintf('Std subsidence (Z): %.3f\n', std(Zc));

%% --- 2. Sampling opsional supaya ringan -------------------------------
Nmax = 2000;  % sesuaikan dengan RAM (2000–3000 cukup)
if numel(Zc) > Nmax
    idx = randperm(numel(Zc), Nmax);
    X = X(idx); Y = Y(idx); Zc = Zc(idx);
end
n = numel(Zc);

%% --- 3. Hitung semua pasangan jarak & kovariansi ----------------------
fprintf('Membentuk pasangan jarak & kovariansi...\n');
pair = nchoosek(1:n,2);
dx = X(pair(:,1)) - X(pair(:,2));
dy = Y(pair(:,1)) - Y(pair(:,2));
dist_km = hypot(dx,dy) / 1000;          % jarak antar titik (km)
cov_ij  = (Zc(pair(:,1))) .* (Zc(pair(:,2)));

fprintf('Total pairs: %.0f\n', numel(cov_ij));

%% --- 4. Plot sebaran mentah seperti slide pertama ---------------------
figure('Color','w');
scatter(dist_km, cov_ij, 2, 'k', 'filled');
xlabel('Distance (km)','FontWeight','bold');
ylabel('Cov [mm^2 yr^{-2}]','FontWeight','bold');
title('Raw pairwise covariance (all data)','FontWeight','bold');
grid on;
ylim([-5000 15000]); xlim([0 100]);
text(40, 12000, 'Subsidence area - Bandung','FontSize',10);

%% --- 5. Windowing untuk menghitung mean tiap jarak --------------------
bin_w = 2;            % window 2 km
max_lag = 100;        % maksimum 100 km
edges   = 0:bin_w:max_lag;
centers = edges(1:end-1) + bin_w/2;
min_pairs = 80;       % minimal pasangan per bin
trim_frac = 0.05;     % trimming 5% ekstrem

ecf = nan(size(centers));
for k = 1:numel(centers)
    inb = dist_km >= edges(k) & dist_km < edges(k+1);
    if sum(inb) >= min_pairs
        v = sort(cov_ij(inb));
        t = floor(trim_frac*numel(v));
        if numel(v) > 2*t, v = v(1+t:end-t); end
        ecf(k) = mean(v,'omitnan');
    end
end

ecf_s = movmean(ecf,3,'omitnan');  % haluskan ringan

%% --- 6. Overlay hasil rata-rata (garis merah) -------------------------
hold on;
plot(centers, ecf_s, 'r-', 'LineWidth',2);
legend({'All pairs','Mean per 2 km window'},'Location','northeast');

%% --- 7. Plot khusus ECF halus (mirip slide terakhir) ------------------
figure('Color','w');
plot(centers, ecf_s, '-o','LineWidth',1.8,'Color',[0 0.447 0.741]);
xlabel('Distance (km)','FontWeight','bold');
ylabel('Cov [mm^2 yr^{-2}]','FontWeight','bold');
title('Empirical Covariance Function of Subsidence (Bandung)','FontWeight','bold');
grid on; yline(0,'k:');

fprintf('\nECF stats -> max: %.2f, mean: %.2f, min: %.2f\n', ...
        max(ecf_s,[],'omitnan'), mean(ecf_s,'omitnan'), min(ecf_s,[],'omitnan'));
fprintf('Bin width : %.0f km | Max lag : %.0f km | Pairs/bin ≥ %d\n', ...
        bin_w, max_lag, min_pairs);

%% ==========================================================
%  TASK C - POINT 3 : Analytical Fit (Hirvonen & Exponential)
%  ==========================================================

fprintf('\n=== STEP 2: Analytical Fit (Hirvonen & Exponential) ===\n');

ok = isfinite(ecf_s);
xfit = centers(ok);
yfit = ecf_s(ok);

% --- Normalisasi agar stabil
yfit = yfit - mean(yfit,'omitnan');
if max(abs(yfit)) > 1e4
    fprintf('Normalizing (dividing by 1e4)...\n');
    yfit = yfit / 1e4;
elseif max(abs(yfit)) > 1e3
    fprintf('Normalizing (dividing by 1e3)...\n');
    yfit = yfit / 1e3;
end

% --- Definisikan fungsi analitik
hirv = @(p,h) p(1)./(1+(h./p(2))).^p(3);   % [C0, d0, p]
expf = @(p,h) p(1)*exp(-h./p(2));          % [C0, a]

% --- Tebakan awal
p0_h = [max(yfit),5,1.2];
p0_e = [max(yfit),5];

% --- Fungsi objektif (MSE)
obj_h = @(p) mean((yfit - hirv(abs(p),xfit)).^2);
obj_e = @(p) mean((yfit - expf(abs(p),xfit)).^2);

% --- Optimasi manual (tanpa toolbox)
ph = fminsearch(obj_h,p0_h,optimset('Display','off','MaxFunEvals',5000));
pe = fminsearch(obj_e,p0_e,optimset('Display','off','MaxFunEvals',5000));

% --- Evaluasi hasil
y_h = hirv(abs(ph),xfit);
y_e = expf(abs(pe),xfit);
rmse_h = sqrt(mean((yfit - y_h).^2));
rmse_e = sqrt(mean((yfit - y_e).^2));

fprintf('HIRVONEN:  C0=%.2f, d0=%.2f km, p=%.2f, RMSE=%.2f\n',abs(ph(1)),abs(ph(2)),abs(ph(3)),rmse_h);
fprintf('EXPONENTIAL:  C0=%.2f, a=%.2f km, RMSE=%.2f\n',abs(pe(1)),abs(pe(2)),rmse_e);

% --- Jarak di mana C ~ 5% C0
thr = 0.05*abs(ph(1));
d_no = fminsearch(@(h) abs(hirv(abs(ph),h)-thr),10);
fprintf('Distance where cov ~5%% of C0 (Hirvonen): %.1f km\n',d_no);

% --- Plot hasil
hh = linspace(0,max(centers),300);
figure('Color','w'); hold on; grid on;
plot(xfit,yfit,'bo-','LineWidth',1.3,'DisplayName','Empirical Covariance (ECF)');
plot(hh,hirv(abs(ph),hh),'r-','LineWidth',2.2,'DisplayName','Hirvonen fit');
plot(hh,expf(abs(pe),hh),'k--','LineWidth',1.5,'DisplayName','Exponential fit');
xlabel('Distance (km)','FontWeight','bold');
ylabel('Cov [mm^2 yr^{-2}]','FontWeight','bold');
title('Empirical vs Analytical Covariance Models','FontWeight','bold');
legend('Location','northeast');
yline(0,'k:');
text(d_no,thr,sprintf('~%.1f km (no correlation)',d_no), ...
     'FontSize',9,'Color',[0.3 0.3 0.3],'VerticalAlignment','bottom');


%% ============================================================
%  TASK C - Point 4 : Ordinary Kriging (Adaptive & Stable)
%  ============================================================
clear; clc; close all;

fprintf('\n=== STEP 3: Ordinary Kriging Surface and Error ===\n');

% ------------------------------------------------------------
% 1. Load data (Bandung Subsidence)
% ------------------------------------------------------------
load('data.mat');

% Pastikan orientasi kolom (3039×1)
E = double(data.east(:));
N = double(data.north(:));
Z = double(data.DV(:));

% Hapus NaN (kalau ada)
ok = isfinite(E) & isfinite(N) & isfinite(Z);
E = E(ok); N = N(ok); Z = Z(ok);

% Gabung jadi koordinat XY
XY = [E N];
fprintf('Jumlah data: %d titik\n', numel(Z));

% ------------------------------------------------------------
% 2. AOI (Area of Interest) dan Grid 250 m
% ------------------------------------------------------------
x_min = 770000; x_max = 813000;
y_min = 9215000; y_max = 9248000;
[xq, yq] = meshgrid(x_min:250:x_max, y_min:250:y_max);

fprintf('Grid size: %d × %d nodes\n', size(xq,1), size(xq,2));

% ------------------------------------------------------------
% 3. Variogram Model (Exponential)
% ------------------------------------------------------------
sigmaZ = std(Z,'omitnan');
C0  = sigmaZ^2;         % sill (variance)
a   = 8000;             % range ≈ 8 km (aktif)
nug = 0.05 * C0;        % nugget 5% sill

% Semivariogram model
gamma = @(h) nug + C0*(1 - exp(-h./a));
fprintf('Variogram params:  C0=%.2f | a=%.0f m | nug=%.2f\n', C0, a, nug);

% ------------------------------------------------------------
% 4. Hapus titik duplikat
% ------------------------------------------------------------
[XYu, ia] = unique(XY, 'rows');
Zu = Z(ia);
nu = numel(Zu);
fprintf('Data unik: %d titik\n', nu);

% ------------------------------------------------------------
% 5. Ordinary Kriging (Local Adaptive)
% ------------------------------------------------------------
Kmin = 12;       % minimal titik tetangga
Kmax = 60;       % maksimal titik tetangga
R0   = 3000;     % radius awal 3 km
Rmax = 30000;    % radius maksimum 30 km
epsR = 1e-6 * C0; % regularisasi kecil (stabilitas numerik)

Zgrid   = nan(size(xq));
SEEgrid = nan(size(xq));

fprintf('Running adaptive local ordinary kriging ...\n');

for i = 1:numel(xq)
    % Hitung jarak ke semua titik data
    dAll = hypot(XYu(:,1) - xq(i), XYu(:,2) - yq(i));
    I = find(dAll <= R0);
    R = R0;
    
    % Perbesar radius hingga dapat minimal Kmin titik
    while numel(I) < Kmin && R < Rmax
        R = R * 1.5;
        I = find(dAll <= R);
    end
    if isempty(I)
        continue
    end
    if numel(I) > Kmax
        [~,o] = sort(dAll(I));
        I = I(o(1:Kmax));
    end

    k = numel(I);
    Xl = XYu(I,:);  Zl = Zu(I);
    dl = dAll(I);

    % Matriks semivariogram antar titik
    Dx = Xl(:,1); Dy = Xl(:,2);
    Dl = sqrt((Dx - Dx').^2 + (Dy - Dy').^2);
    G  = gamma(Dl);
    G(1:k+1:end) = epsR;        % stabilisasi diagonal
    G  = [G, ones(k,1); ones(1,k), 0];

    % Vektor semivariogram ke grid node
    g = gamma(dl(:));
    g = [g; 1];

    % Bobot kriging (pakai pseudo-inverse)
    w = pinv(G) * g;
    lam = w(1:k);
    mu  = w(end);

    % Estimasi dan standard error
    Zgrid(i)   = lam' * Zl;
    SEEgrid(i) = sqrt(max(lam' * g(1:k) + mu, 0));
end

% ------------------------------------------------------------
% 6. Statistik Hasil
% ------------------------------------------------------------
fprintf('\n=== STATISTICS ===\n');
fprintf('Scattered data: mean=%6.2f | std=%6.2f | min=%6.2f | max=%6.2f\n', ...
        mean(Z,'omitnan'), std(Z,'omitnan'), min(Z), max(Z));
fprintf('Kriged surface: mean=%6.2f | std=%6.2f | min=%6.2f | max=%6.2f\n', ...
        mean(Zgrid(:),'omitnan'), std(Zgrid(:),'omitnan'), ...
        min(Zgrid(:)), max(Zgrid(:)));
fprintf('SEEgrid:        mean=%6.2f | std=%6.2f | min=%6.2f | max=%6.2f\n', ...
        mean(SEEgrid(:),'omitnan'), std(SEEgrid(:),'omitnan'), ...
        min(SEEgrid(:)), max(SEEgrid(:)));

% ------------------------------------------------------------
% 7. Plot Hasil
% ------------------------------------------------------------
figure('Color','w');
imagesc(x_min:x_max, y_min:y_max, Zgrid);
set(gca,'YDir','normal'); axis equal tight;
colormap(parula); colorbar;
caxis([min(Z) max(Z)]);
xlabel('Easting (m)'); ylabel('Northing (m)');
title('Ordinary Kriging Subsidence Surface (Bandung)');

figure('Color','w');
imagesc(x_min:x_max, y_min:y_max, SEEgrid);
set(gca,'YDir','normal'); axis equal tight;
colormap(turbo); colorbar;
xlabel('Easting (m)'); ylabel('Northing (m)');
title('Standard Error of Estimates (Kriging)');
