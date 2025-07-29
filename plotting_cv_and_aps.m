clear all
close all
clc

% Define the r_pn values (in Ohm·cm or whatever unit you're using)
rpn_grid = [50, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 2200, 2400, 2600, 2800, 3000];
Cmy_grid = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1, 0.11, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.19, 0.2];
Rmy_grid = [50, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000];
Ri_grid = [0.010, 0.020, 0.030, 0.040, 0.050, 0.060, 0.070, 0.080, 0.090, 0.100, 0.110, 0.120, 0.130, 0.140, 0.150, 0.160, 0.170, 0.180, 0.190, 0.200];
HH_Temp_grid = [6.3, 7, 9, 11, 13, 15, 17, 19, 21, 23, 25, 27, 29, 31];
SC_Temp_grid = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45];
DC_Temp_grid = [20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45];


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

Ri_cv_set2_w_set1rpn = [0 0 0 0 0 0	0	0.0996121419035100	0.0952646417217168	0.0901850315354267	0.0856116806717296	0.0815736512994839	0.0782062864631637	0.0751300045584766	0.0725338102572954	0.0700579074006869	0.0680573413781243	0.0661995143674720	0.0645365519967455	0.0630358058388596];
Ri_ap_duration_set2_w_set1rpn = [0 0 0 0 0 0 0	17.4800000000000	16.4000000000000	16.1000000000000	15.9600000000000	15.9000000000000	15.8800000000000	15.8600000000000	15.8600000000000	15.8400000000000	15.8200000000000	15.8200000000000	15.7800000000000	15.7600000000000];
Ri_ap_peak_set2_w_set1rpn = [-56.0750080663227	-57.2915524240940	-57.4496713097949	-57.5406212166057	-57.5754554649410	-57.5207469331020	-56.8230369905049	-28.7654694015246	-26.5758562528112	-26.0341746062314	-25.8298823860426	-25.7445994878042	-25.7092971558443	-25.6955075056898	-25.6899008128874	-25.6854301862295	-25.6783614008883	-25.6665585992422	-25.6488539461051	-25.6245709238854];
Ri_cv_set1_w_set2Cmy = [0 0 0 0 0 0 0 0 0	0.0983040402543857	0.0893354704798536	0.0842288005290468	0.0808167920083946	0.0781421177348171	0.0759072221808596	0.0739565472908702	0.0721320756776698	0.0705044057241812	0.0689908452203534	0.0675645781575313];
Ri_ap_duration_set1_w_set2Cmy = [0 0 0 0 0 0 0 0 0	18.2400000000000	17.2000000000000	17.1000000000000	17.0800000000000	17.0400000000000	17.0400000000000	17	16.9800000000000	17	16.9800000000000	17.0200000000000];
Ri_ap_peak_set1_w_set2Cmy = [-56.4700302942053	-57.4253447716838	-57.5859476211292	-57.6914129915418	-57.7616153760346	-57.8034661824813	-57.8136719247623	-57.7624774746022	-57.3258749369026	-50.5615355466535	-44.6979870222734	-42.7351058926109	-42.1710973054186	-42.1756349636920	-42.4568605870496	-42.8887016065210	-43.4069057651976	-43.9743785106935	-44.5675596392011	-45.1704243051615];

% T_actual Data
HH_Temp_cv = [13.6073573573574	13.9880952380952	15.0401069518717	16.1290322580645	17.0977011494253	18.3531746031746	19.4230769230769	20.6250000000000	21.7391304347826	22.4802371541502	23.5389610389610	23.8095238095238	24.1071428571429	22.9978354978355];
HH_Temp_ap_duration = [3.20000000000000	2.98000000000000	2.44000000000000	2	1.66000000000000	1.38000000000000	1.16000000000000	0.980000000000000	0.840000000000000	0.740000000000000	0.640000000000000	0.580000000000000	0.540000000000000	0.520000000000000];
HH_Temp_ap_peak = [37.5644082169696	37.1362116187897	35.7457022475952	34.0573081686096	32.0407279761778	29.6329070576377	26.7647919493737	23.4370603278342	19.5287422773255	14.9830011318067	9.71508116863436	3.67226033646928	-3.32523772583959	-11.9224614482906];
SC_Temp_cv = [0.744708994708995	0.772272024540627	0.799858231437179	0.826335470085470	0.860784313725491	0.867582875229053	0.888687600644123	0.933848286789463	0.959353146853147	0.970308123249300	1.00892857142857	1.01150895140665	1.03782474088033	1.06608851674641	1.09153091060986	1.12249066293184	1.15721288515406	1.16911764705882	1.18227554179567	1.20065789473684	1.24232456140351	1.25694444444444	1.28075396825397 1.30158730158730	1.30456349206349	1.31792717086835];
SC_Temp_ap_duration = [4.64000000000000	4.30000000000000	3.98000000000000	3.70000000000000	3.44000000000000	3.20000000000000	2.98000000000000	2.78000000000000	2.60000000000000	2.42000000000000	2.26000000000000	2.12000000000000	1.98000000000000	1.88000000000000	1.74000000000000	1.66000000000000	1.56000000000000	1.48000000000000	1.40000000000000	1.32000000000000	1.26000000000000	1.22000000000000	1.16000000000000 1.12000000000000	1.08000000000000	1.04000000000000];
SC_Temp_ap_peak = [23.4455921786601	23.0182860890937	22.5600064154362	22.0669569492862	21.5400066660028	20.9735428720897	20.3674031778080	19.7217998467669	19.0297677110693	18.2904831265601	17.5039876839092	16.6659650648531	15.7727385967574	14.8285525400549	13.8247617461328	12.7624940974415	11.6300331988821	10.4516357396269	9.20204308735702	7.88447640431262	6.49207634086779	5.04567671119379	3.50853136002335 1.92164101608220	0.252800082603176	-1.49752152103146];
DC_Temp_cv = [0.248688529938530	0.254809616906391	0.260450217641998	0.266184348376129	0.271476318872679	0.276245280720984	0.279046401169163	0.283024106963411	0.282661298957151	0.282661298957151	0.279108658388729	0.269435085968004	0.252872883456118	0.213171612081267, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
DC_Temp_ap_duration = [7.30000000000000	6.98000000000000	6.70000000000000	6.44000000000000	6.22000000000000	6.02000000000000	5.86000000000000	5.74000000000000	5.64000000000000	5.62000000000000	5.64000000000000	5.76000000000000	6.02000000000000	6.40000000000000, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0];
DC_Temp_ap_peak = [1.98697024137971	0.417039747088619	-1.21431020737798	-2.91018180649520	-4.67424648983438	-6.51255088888750	-8.43374838643677	-10.4507298485085	-12.5833199035873	-14.8626292030437	-17.3392571687887	-20.1022245176226	-23.3304197256292	-27.4349764022747, -33.5212089083390	-42.6519345852629	-49.7997625398069	-53.2945942373314	-55.0058116412702	-55.9355136906158	-56.4936715134550	-56.8565550225845	-57.1072137187832	-57.2885807111151	-57.4245858471018	-57.5294355440625];

SC_Cmy_Temp33_cv = [1.67249417249417	1.67249417249417	1.60256410256410	1.46825396825397	1.38437001594896	1.34270334928230	1.25961538461539	1.21094591682827	1.16713352007470	1.13250148544266	1.08070839978735	1.03782474088033	1.02148487159929	0.995169082125604	0.956739675275144	0.946428571428571	0.922619047619048	0.903463203463204	0.892640692640693	0.874125874125874];
SC_Cmy_Temp33_ap_duration = [1.86000000000000	1.86000000000000	1.86000000000000	1.88000000000000	1.86000000000000	1.86000000000000	1.86000000000000	1.88000000000000	1.88000000000000	1.86000000000000	1.86000000000000	1.86000000000000	1.88000000000000	1.88000000000000	1.86000000000000	1.86000000000000	1.86000000000000	1.86000000000000	1.86000000000000	1.86000000000000];
SC_Cmy_Temp33_ap_peak = [9.35611938494618	12.4810495785993	13.4869672165209	13.9698265005254	14.2559141313523	14.4329224000717	14.5639356426521	14.6545913704830	14.7208273776724	14.7721817883016	14.8140509433859	14.8490130754848	14.8776867549625	14.8993285098819	14.9122687074704	14.9322949480972	14.9439146676756	14.9542805831562	14.9630801156516	14.9718869656128];
DC_Cmy_Temp33_cv = [0 0 0 0	0.110985837346441	0.154178326145246	0.167535068525168	0.187028967369044	0.208985044038663	0.218000758659050	0.233478964712227	0.257500384496245	0.296980796980797	0.324487525822342	0.346631656976485	0.365747314039534	0.383167018006435	0.396308159321960	0.406339610943192	0.416855631141346];
DC_Cmy_Temp33_ap_duration = [0 0 0 0	15.1600000000000	13.4800000000000	11.6400000000000	9.66000000000000	7.70000000000000	6.62000000000000	6.12000000000000	5.72000000000000	4.98000000000000	4.48000000000000	4.12000000000000	3.82000000000000	3.62000000000000	3.44000000000000	3.28000000000000	3.16000000000000];
DC_Cmy_Temp33_ap_peak = [-58.1124000000000	-58.1124000000000	-57.9321139947273	-57.5726101708331	-57.0732866129074	-56.3859784792385	-55.2940526916373	-53.1812675804103	-48.2501143819353	-38.0250507974745	-29.2051997120745	-24.1546801012557	-20.6566349458784	-17.9326257003798	-15.6818163066244	-13.7626640398896	-12.0881179613967	-10.6071059984017	-9.28351271130768	-8.09034555415389];


% % PLOTTING r_pn
% %%%%%%%%%%%%%%%
% figure(1);
% plot(rpn_grid, rpn_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(rpn_grid, rpn_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(2);
% plot(rpn_grid, rpn_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(rpn_grid, rpn_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(3);
% plot(rpn_grid, rpn_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(rpn_grid, rpn_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$r_{pn}$ ($10^6$ k$\Omega$/cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $R_i = 0.0712$', 'Set (2) $C_{my} = 0.0379$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% % PLOTTING C_my
% %%%%%%%%%%%%%%%
% figure(4);
% plot(Cmy_grid, Cmy_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, Cmy_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(5);
% plot(Cmy_grid, Cmy_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, Cmy_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(6);
% plot(Cmy_grid, Cmy_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, Cmy_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $r_{pn} = 321\cdot10^6$, $R_i = 0.0712$', 'Set (2) $r_{pn} = 2450\cdot10^6$, $R_i = 0.155$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% % PLOTTING R_my
% %%%%%%%%%%%%%%%
% figure(7);
% plot(Rmy_grid, Rmy_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Rmy_grid, Rmy_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(8);
% plot(Rmy_grid, Rmy_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Rmy_grid, Rmy_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(9);
% plot(Rmy_grid, Rmy_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Rmy_grid, Rmy_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{my}$ (k$\Omega\cdot$cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1)', 'Set (2)', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% % PLOTTING R_i
% %%%%%%%%%%%%%%
% figure(10);
% plot(Ri_grid, Ri_cv_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Ri_grid, Ri_cv_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(11);
% plot(Ri_grid, Ri_ap_peak_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Ri_grid, Ri_ap_peak_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% figure(12);
% plot(Ri_grid, Ri_ap_duration_set1, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on
% plot(Ri_grid, Ri_ap_duration_set2, '-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
% legend('Set (1) w/ $C_{my} = 0.113$, $r_{pn} = 321\cdot10^6$', 'Set (2) $C_{my} = 0.0379$, $r_{pn} = 2450\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% % grid on;
% 
% PLOTTING R_i (w/ Set (2) params and r_pn = 321*10^6)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(13);
plot(Ri_grid, Ri_cv_set1_w_set2Cmy, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on
plot(Ri_grid, Ri_cv_set2_w_set1rpn, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) w/ $C_{my} = 0.0379$', 'Set (2) w/ $r_{pn} = 321\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(14);
plot(Ri_grid, Ri_ap_peak_set1_w_set2Cmy, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Ri_grid, Ri_ap_peak_set2_w_set1rpn, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) w/ $C_{my} = 0.0379$', 'Set (2) w/ $r_{pn} = 321\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
set(gca, 'FontSize', 13);
% grid on;

figure(15);
plot(Ri_grid, Ri_ap_duration_set1_w_set2Cmy, '-o', 'LineWidth', 2, 'MarkerSize', 8);
hold on;
plot(Ri_grid, Ri_ap_duration_set2_w_set1rpn, '-o', 'LineWidth', 2, 'MarkerSize', 8);
ylabel('AP Duration (mV)', 'Interpreter', 'latex', 'FontSize', 16)
xlabel('$R_{i}$ (k$\Omega$ cm)', 'Interpreter', 'latex', 'FontSize', 16)
legend('Set (1) w/ $C_{my} = 0.0379$', 'Set (2) w/ $r_{pn} = 321\cdot10^6$', 'Interpreter', 'latex', 'FontSize', 16)
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

% % SC and DC Cmy Analysis for T_actual = 33
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure(22);
% plot(Cmy_grid, SC_Cmy_Temp33_cv, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, DC_Cmy_Temp33_cv, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Conduction Velocity (m/s)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% title('CV vs C_{my}')
% legend('SC: $T_{actual} = 33^\circ$C', 'DC: $T_{actual} = 33^\circ$C', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(23);
% plot(Cmy_grid, SC_Cmy_Temp33_ap_duration, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, DC_Cmy_Temp33_ap_duration, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('AP Duration (ms)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% title('AP Duration vs C_{my}')
% legend('SC: $T_{actual} = 33^\circ$C', 'DC: $T_{actual} = 33^\circ$C', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;
% 
% figure(24);
% plot(Cmy_grid, SC_Cmy_Temp33_ap_peak, 'r-o', 'LineWidth', 2, 'MarkerSize', 8);
% hold on;
% plot(Cmy_grid, DC_Cmy_Temp33_ap_peak, 'g-o', 'LineWidth', 2, 'MarkerSize', 8);
% ylabel('Peak Voltage of AP (mV)', 'Interpreter', 'latex', 'FontSize', 16)
% xlabel('$C_{my}$ ($\mu$F/cm$^2$)', 'Interpreter', 'latex', 'FontSize', 16)
% title('Peak Voltage of AP vs C_{my}')
% legend('SC: $T_{actual} = 33^\circ$C', 'DC: $T_{actual} = 33^\circ$C', 'Interpreter', 'latex', 'FontSize', 16)
% set(gca, 'FontSize', 13);
% grid on;

