clear all
close all
clc

% Define the r_pn values (in Ohm·cm or whatever unit you're using)
rpn_grid = [50, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000];
Cmy_grid = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];
Rmy_grid = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
Ri_grid = [0.020, 0.040, 0.060, 0.080, 0.100, 0.120, 0.140, 0.160, 0.180, 0.200];
HH_Temp_grid = [6.3, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31];
SC_Temp_grid = [20, 22, 24, 26, 28, 30, 32, 34, 36, 38, 40, 42, 44, 46, 48, 50, 52, 54, 56];
DC_Temp_grid = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33];


% r_pn Data
rpn_cv_set1 = [0.227832272223429	0.238836009669343	0.254338563192245	0.270114864215988	0.284683082071032	0.299832134736560	0.312945156695157	0.329279041071743	0.341255868544601	0.353573213393303	0.366834657879434	0.379472379472380	0.388755464795214	0.403043351335571	0.410990712074304	0.421921714173530];
rpn_ap_duration_set1 = [7.40000000000000	7.34000000000000	7.24000000000000	7.16000000000000	7.02000000000000	6.88000000000000	6.76000000000000	6.62000000000000	6.48000000000000	6.36000000000000	6.24000000000000	6.12000000000000	6	5.90000000000000	5.80000000000000	5.70000000000000];
rpn_ap_peak_set1 = [1.41966294247315	1.71766407439354	2.17571160982705	2.69241327074242	3.25347915996051	3.84466780624659	4.45205823082682	5.06367674427982	5.66999084632803	6.26470493545118	6.84280662103053	7.40039675138128	7.93542070508890	8.44729613004956	8.93558815104667	9.40029193431265];
rpn_cv_set2 = [0.0688049032222233	0.0721256920575004	0.0767038094997734	0.0817234845602019	0.0868664591778791	0.0923529408395039	0.0978094342368883	0.103584764146517	0.109743538959151	0.115918187431086	0.122227730498407	0.128570545237212	0.134954762094797	0.141532518565909	0.147743173484630	0.154451205508560];
rpn_ap_duration_set2 = [16.0200000000000	15.7000000000000	15.2800000000000	14.8400000000000	14.3800000000000	13.9400000000000	13.4800000000000	13.0200000000000	12.5800000000000	12.1200000000000	11.6800000000000	11.2400000000000	10.8200000000000	10.4200000000000	10.0200000000000	9.64000000000000];
rpn_ap_peak_set2 = [-25.5706110991875	-24.9309985293710	-24.0511683392359	-23.1475825830547	-22.2263022109718	-21.2921286850994	-20.3489408289470	-19.4002809167801	-18.4497262037669	-17.5005941853647	-16.5560900237029	-15.6196197918003	-14.6940724712363	-13.7825869281457	-12.8878289041651	-12.0114531696413];

% C_my Data
Cmy_cv_set1 = [0 0 0 0 0.138972714445893	0.168824276030921	0.190797678682091	0.208114028092337	0.222982408059927	0.234939434710022	0.245606761762314	0.254809616906391	0.261113225901958	0.268823921372647	0.274526089505570	0.280143855076099	0.283792009043085	0.289035913414023	0.293305652680653	0.295262293086660];
Cmy_ap_duration_set1 = [0 0 0 0 13.4200000000000	10.9000000000000	9.58000000000000	8.76000000000000	8.18000000000000	7.72000000000000	7.38000000000000	7.12000000000000	6.88000000000000	6.70000000000000	6.54000000000000	6.40000000000000	6.30000000000000	6.18000000000000	6.08000000000000	6];
Cmy_ap_peak_set1 = [-58.1124000000000	-58.1124000000000	-58.0438032053246	-57.6731333605780	-19.7725438734397	-13.0256203796846	-8.60259984706558	-5.24721782101643	-2.56820572245506	-0.364220050688241	1.48735922173004	3.06804104535706	4.43520281173143	5.63032274620713	6.68514021222542	7.62279433803886	8.46330887284787	9.22016155083013	9.90570255671926	10.5302987552505];
Cmy_cv_set2 = [0, 0, 0.0767571713086678	0.146075465136498	0.181719227062257	0.206219199684605	0.224238551633497	0.237891054841427	0.248072383949646	0.255666076650762	0.261440808126092	0.266615737203973	0.270082815734990	0.272900810402736	0.274701099952176	0.276354786158708	0.276551339164743	0.277452240065644	0.277648793071679	0.277322944103014];
Cmy_ap_duration_set2 = [0, 0, 16.3600000000000	10.0600000000000	8.22000000000000	7.28000000000000	6.74000000000000	6.36000000000000	6.10000000000000	5.90000000000000	5.76000000000000	5.64000000000000	5.52000000000000	5.44000000000000	5.36000000000000	5.32000000000000	5.24000000000000	5.20000000000000	5.16000000000000	5.12000000000000];
Cmy_ap_peak_set2 = [-58.3540000000000	-57.8149440381218	-28.4546715076982	-12.4325542943729	-5.42863059221648	-0.795186679454457	2.57627833416860	5.15833607003834	7.20509797915812	8.86940835051960	10.2499954334081	11.4143494481623	12.4091935118673	13.2694209495534	14.0201025233308	14.6798581581332	15.2659106364325	15.7887226638238	16.2579690161492	16.6807621924412];

% R_my Data
Rmy_cv_set1 = [0.247877105076472	0.248688529938530	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507	0.248425350690507];
Rmy_ap_duration_set1 = [7.30000000000000	7.30000000000000	7.30000000000000	7.30000000000000	7.28000000000000	7.28000000000000	7.28000000000000	7.28000000000000	7.28000000000000	7.28000000000000	7.28000000000000];
Rmy_ap_peak_set1 = [1.98325356126044	1.99178694275012	1.99591176729788	1.99726551774156	1.99793840489875	1.99834085911649	1.99860862882319	1.99879963189528	1.99894274126072	1.99905396380672	1.99914288849890];
Rmy_cv_set2 = [0.134308395006141	0.135882189132963	0.136373583211761	0.136822589740622	0.136755632763592	0.136822589740622	0.137040435679974	0.137040435679974	0.137040435679974	0.137103141867207	0.136818338950825];
Rmy_ap_duration_set2 = [10.8600000000000	10.7800000000000	10.7200000000000	10.7200000000000	10.7200000000000	10.7000000000000	10.7000000000000	10.7000000000000	10.7000000000000	10.7000000000000	10.6800000000000];
Rmy_ap_peak_set2 = [-14.7679914844735	-14.5778379164774	-14.4811065818496	-14.4486133091946	-14.4322451393654	-14.4226105139786	-14.4159942317130	-14.4112590212180	-14.4077285372322	-14.4049950039345	-14.4028159909922];

% R_i Data
Ri_cv_set1 = [0.656273620559335	0.365321858753978	0.276379029004704	0.231697133862666	0.203865210118467	0.185144264602393	0.170439159678547	0.159229465814205	0.150116410009450	0.142184205446954];
Ri_ap_duration_set1 = [7.26000000000000	7.32000000000000	7.30000000000000	7.28000000000000	7.26000000000000	7.28000000000000	7.26000000000000	7.24000000000000	7.22000000000000	7.22000000000000];
Ri_ap_peak_set1 = [2.15660218526377	2.02788415504866	1.99425656575230	1.98558415011015	1.99434054351257	2.01529677045871	2.04373934931826	2.07583588780410	2.10858147671473	2.14038971824492];
Ri_cv_set2 = [0	0.397064579256359	0.266667547096208	0.207453978966915	0.176085390656499	0.157485041922638	0.144221473611932	0.134378008340077	0.126154036566979	0.119651325068391];
Ri_ap_duration_set2 = [0	11.3000000000000	10.9200000000000	10.8800000000000	10.8400000000000	10.7800000000000	10.7600000000000	10.7200000000000	10.6800000000000	10.6400000000000];
Ri_ap_peak_set2 = [-56.0395221479148	-16.2419190112514	-14.5350504219759	-14.4730032914566	-14.4852595753916	-14.4915854356809	-14.4818982187532	-14.4576422635705	-14.4214572562469	-14.3764806198489];

% T_actual Data
HH_Temp_cv = [13.6073573573574	13.9880952380952	15.0401069518717	16.1290322580645	17.0977011494253	18.3531746031746	19.4230769230769	20.6250000000000	21.7391304347826	22.4802371541502	23.5389610389610	23.8095238095238	24.1071428571429	22.9978354978355];
HH_Temp_ap_duration = [3.20000000000000	2.98000000000000	2.44000000000000	2	1.66000000000000	1.38000000000000	1.16000000000000	0.980000000000000	0.840000000000000	0.740000000000000	0.640000000000000	0.580000000000000	0.540000000000000	0.520000000000000];
HH_Temp_ap_peak = [37.5644082169696	37.1362116187897	35.7457022475952	34.0573081686096	32.0407279761778	29.6329070576377	26.7647919493737	23.4370603278342	19.5287422773255	14.9830011318067	9.71508116863436	3.67226033646928	-3.32523772583959	-11.9224614482906];
SC_Temp_cv = [0.744708994708995	0.799858231437179	0.860784313725491	0.888687600644123	0.959353146853147	1.00892857142857	1.03782474088033	1.09153091060986	1.15721288515406	1.18227554179567	1.24232456140351	1.28075396825397	1.30456349206349	1.34173669467787	1.33928571428571	1.34837588881707	1.32090336134454	1.23856209150327	1.08838383838384];
SC_Temp_ap_duration = [4.64000000000000	3.98000000000000	3.44000000000000	2.98000000000000	2.60000000000000	2.26000000000000	1.98000000000000	1.74000000000000	1.56000000000000	1.40000000000000	1.26000000000000	1.16000000000000	1.08000000000000	1.02000000000000	0.960000000000000	0.940000000000000	0.940000000000000	0.980000000000000	1.08000000000000];
SC_Temp_ap_peak = [23.4455921786601	22.5600064154362	21.5400066660028	20.3674031778080	19.0297677110693	17.5039876839092	15.7727385967574	13.8247617461328	11.6300331988821	9.20204308735702	6.49207634086779	3.50853136002335	0.252800082603176	-3.32626634769536	-7.21715359290324	-11.4504622712077	-16.1258221521560	-21.3597584556493	-27.5844883147625];
DC_Temp_cv = [0.248688529938530	0.254809616906391	0.260450217641998	0.266184348376129	0.271476318872679	0.276245280720984	0.279046401169163	0.283024106963411	0.282661298957151	0.282661298957151	0.279108658388729	0.269435085968004	0.252872883456118	0.213171612081267];
DC_Temp_ap_duration = [7.30000000000000	6.98000000000000	6.70000000000000	6.44000000000000	6.22000000000000	6.02000000000000	5.86000000000000	5.74000000000000	5.64000000000000	5.62000000000000	5.64000000000000	5.76000000000000	6.02000000000000	6.40000000000000];
DC_Temp_ap_peak = [1.98697024137971	0.417039747088619	-1.21431020737798	-2.91018180649520	-4.67424648983438	-6.51255088888750	-8.43374838643677	-10.4507298485085	-12.5833199035873	-14.8626292030437	-17.3392571687887	-20.1022245176226	-23.3304197256292	-27.4349764022747];




% PLOTTING r_pn
%%%%%%%%%%%%%%%
figure(1);
plot(rpn_grid, rpn_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(rpn_grid, rpn_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(2);
plot(rpn_grid, rpn_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(rpn_grid, rpn_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(3);
plot(rpn_grid, rpn_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(rpn_grid, rpn_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

% PLOTTING C_my
%%%%%%%%%%%%%%%
figure(4);
plot(Cmy_grid, Cmy_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Cmy_grid, Cmy_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(5);
plot(Cmy_grid, Cmy_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Cmy_grid, Cmy_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(6);
plot(Cmy_grid, Cmy_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Cmy_grid, Cmy_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

% PLOTTING R_my
%%%%%%%%%%%%%%%
figure(7);
plot(Rmy_grid, Rmy_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Rmy_grid, Rmy_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(8);
plot(Rmy_grid, Rmy_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Rmy_grid, Rmy_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(9);
plot(Rmy_grid, Rmy_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Rmy_grid, Rmy_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

% PLOTTING R_i
%%%%%%%%%%%%%%
figure(10);
plot(Ri_grid, Ri_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Ri_grid, Ri_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(11);
plot(Ri_grid, Ri_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Ri_grid, Ri_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(12);
plot(Ri_grid, Ri_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Ri_grid, Ri_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

% % PLOTTING TEMPERATURE DATA
% %%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(13);
% plot(HH_Temp_grid, HH_Temp_cv, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('HH: CV vs Temp.')
% % legend('HH Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(14);
% plot(HH_Temp_grid, HH_Temp_ap_peak, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('HH: Peak Votage of AP vs Temp.')
% % legend('HH Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(15);
% plot(HH_Temp_grid, HH_Temp_ap_duration, 'b-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('HH: AP Duration vs Temp.')
% % legend('HH Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(16);
% plot(SC_Temp_grid, SC_Temp_cv, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('SC: CV vs Temp.')
% % legend('SC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(17);
% plot(SC_Temp_grid, SC_Temp_ap_peak, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('SC: Peak Votage of AP vs Temp.')
% % legend('SC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(18);
% plot(SC_Temp_grid, SC_Temp_ap_duration, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('SC: AP Duration vs Temp.')
% % legend('SC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(19);
% plot(DC_Temp_grid, DC_Temp_cv, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('DC: CV vs Temp.')
% % legend('DC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(20);
% plot(DC_Temp_grid, DC_Temp_ap_peak, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('DC: Peak Votage of AP vs Temp.')
% % legend('DC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(21);
% plot(DC_Temp_grid, DC_Temp_ap_duration, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('Temperature ($^\circ$C)', 'Interpreter', 'latex', 'FontSize', 16)
% title('DC: AP Duration vs Temp.')
% % legend('DC Model', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;

