function [subs] = create_anon_subject_structure()
%% modified version of jimmy's code without subject identifiers
subs = struct();

subs(1).subj_id = 's_01';
subs(1).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_mmsMQ0484_orig.img';
subs(1).mpm_t1 = 's2015-06-18_15-22-154522-00001-00224-1_T1w.nii';
subs(1).nas=[9.9898 142.5147 -8.1787];
subs(1).lpa=[-55.7659 49.4636 -26.2089];
subs(1).rpa=[88.7153 62.4787 -29.1394];

subs(1).coords_lhm=[-12.79 -4.778 52.67];
subs(1).coords_lhv=[-6.143 -55.92 -33.65];
subs(1).coords_rhv=[52.07 -53.15 -39.04];
% subs(1).camera_up_vector=dict();
% subs(1).camera_up_vector('topdown')=[-0.06975 0.9974 -0.01745];
% subs(1).camera_up_vector('back')=[0.0,0.7,1];
% subs(1).camera_position=dict();
% subs(1).camera_position('topdown')=[18.4 57.27 1745];
% subs(1).camera_position('back')=[103.152 -1552.014 -736.475];
% 
subs(2).subj_id = 's_02';
subs(2).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-07_14-18-142013-00001-00192-1.nii';
subs(2).mpm_t1 = 's2014-10-07_14-18-145208-00001-00224-1_T1w.nii';
subs(2).nas=[0.4000 117.3122 -9.6584];
subs(2).lpa=[-77.9750 42.3162 -26.4274];
subs(2).rpa=[71.2590 40.7802 -30.1864];

subs(2).coords_lhm=[-41.13 -6.51 59.56];
subs(2).coords_lhv=[-30.45 -71.09 3.36];
subs(2).coords_rhv=[23.64 -67.55 2.84];
% subs(2).camera_up_vector=dict();
% subs(2).camera_up_vector('topdown')=[0.034 0.981 0.191];
% subs(2).camera_up_vector('back')=[0.0,0.0,1];
% subs(2).camera_position=dict();
% subs(2).camera_position('topdown')=[-15.84 -315.94 1685.15];
% subs(2).camera_position('back')=[52.76 -1680.78 -150.43];

subs(3).subj_id = 's_03';
subs(3).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-02_12-11-121418-00001-00192-1.nii';
subs(3).mpm_t1 = 's2014-10-02_12-11-124515-00001-00224-1_T1w.nii';
subs(3).nas=[2.5710 111.5040 20.5150];
subs(3).lpa=[-78.2450 46.0940 -13.5810];
subs(3).rpa=[73.3740 41.3980 -10.9450];

subs(3).coords_lhm=[-43.78 -14.11 76.92];
subs(3).coords_lhv=[-43.06 -70.82 3.352];
subs(3).coords_rhv=[19.78 -73.76 7.572];
% subs(3).camera_up_vector=dict();
% subs(3).camera_up_vector('topdown')=[0.06759 0.975 0.2115];
% subs(3).camera_up_vector('back')=[0.0,0.0,1];
% subs(3).camera_position=dict();
% subs(3).camera_position('topdown')=[-30.93 -352.1 1707];
% subs(3).camera_position('back')=[-109.70 -1645.32 -325.10];

subs(4).subj_id = 's_04';
subs(4).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-08_11-09-111650-00001-00192-1.nii';
subs(4).mpm_t1 = 's2014-10-08_11-09-114810-00001-00224-1_T1w.nii';
subs(4).nas=[-6.7000 120.9000 -14.8000];
subs(4).lpa=[-86 39.4000 -57.6000];
subs(4).rpa=[75.2000 52.2000 -57.4000];

subs(4).coords_lhm=[-40.05 -23.79 24.16];
subs(4).coords_lhv=[-35.11 -63.54 -53.7];
subs(4).coords_rhv=[27.17 -58.49 -55.37];
% subs(4).camera_up_vector=dict();
% subs(4).camera_up_vector('topdown')=[-0.01112 0.8755 0.4831];
% subs(4).camera_up_vector('back')=[0.0,0.0,1];
% subs(4).camera_position=dict();
% subs(4).camera_position('topdown')=[7.056 -843 1525];
% subs(4).camera_position('back')=[32.37 -1532.98 -304.14];

subs(5).subj_id = 's_05';
subs(5).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-10_16-08-161009-00001-00192-1.nii';
subs(5).mpm_t1 = 's2014-10-10_16-08-163926-00001-00224-1_T1w.nii';
subs(5).nas=[-5.5930 109.1525 16.5860];
subs(5).lpa=[-78.5530 33.6745 -11.0340];
subs(5).rpa=[73.7700 38.6545 -10.3130];

subs(5).coords_lhm=[-42.19 -16.16 75.05];
subs(5).coords_lhv=[-25.99 -73.9 -7.039];
subs(5).coords_rhv=[31.14 -74.6 -9.439];
% subs(5).camera_up_vector=dict();
% subs(5).camera_up_vector('topdown')=[-0.04395 0.9878 0.1496];
% subs(5).camera_up_vector('back')=[0.0,-0.5,1];
% subs(5).camera_position=dict();
% subs(5).camera_position('topdown')=[10.36 -252.4 1729];
% subs(5).camera_position('back')=[137.146 -1631.08 146.89];

subs(6).subj_id = 's_06';
subs(6).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-09_14-47-144933-00001-00192-1.nii';
subs(6).mpm_t1 = 's2014-10-09_14-47-151855-00001-00224-1_T1w.nii';
subs(6).nas=[4.6030 109.0934 26.2327];
subs(6).lpa=[-69.8850 36.9034 -3.0143];
subs(6).rpa=[71.8630 31.1924 -2.7183];

subs(6).coords_lhm=[-37.26 -18.04 84.5];
subs(6).coords_lhv=[-30.67 -73.99 13.46];
subs(6).coords_rhv=[23.78 -74.78 15.74];
% subs(6).camera_up_vector=dict();
% subs(6).camera_up_vector('topdown')=[0.04067 0.9351 0.3521];
% subs(6).camera_up_vector('back')=[0.0,0.0,1];
% subs(6).camera_position=dict();
% subs(6).camera_position('topdown')=[-27.01 -598 1639];
% subs(6).camera_position('back')=[-33.77 -1599.02 -193.16];

subs(7).subj_id = 's_07';
subs(7).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2015-06-09_11-44-114620-00001-00192-1.nii';
subs(7).mpm_t1 = 's2015-06-09_11-44-121745-00001-00224-1_T1w.nii';
subs(7).nas=[-9.4100 110.0300 9.3589];
subs(7).lpa=[-74.5140 46.7120 -17.2951];
subs(7).rpa=[61.0680 49.7940 -11.8531];

subs(7).coords_lhm=[-39.07 -14.66 68.16];
subs(7).coords_lhv=[-26.94 -72.86 -2.085];
subs(7).coords_rhv=[21.79 -71.65 -3.438];
% subs(7).camera_up_vector=dict();
% subs(7).camera_up_vector('topdown')=[-0.03274 0.9577 0.2859];
% subs(7).camera_up_vector('back')=[0.0,-0.7,1];
% subs(7).camera_position=dict();
% subs(7).camera_position('topdown')=[11.44 -484.1 1677];
% subs(7).camera_position('back')=[111.15 -1663.22 -202.43];

subs(8).subj_id = 's_08';
subs(8).headcast_t1 = 'pd_nw_mtflash3d_v3e_HeadCast_12ch_0002/anon_s2014-10-03_15-14-151946-00001-00192-1.nii';
subs(8).mpm_t1 = 's2014-10-03_15-14-155029-00001-00224-1_T1w.nii';
subs(8).nas=[12.9130 110.7524 9.6484];
subs(8).lpa=[-71.1920 37.3994 -21.4386];
subs(8).rpa=[81.9700 24.4084 -23.5716];

subs(8).coords_lhm=[-33.06 -18.01 73.6];
subs(8).coords_lhv=[-29.74 -79.02 -6.597];
subs(8).coords_rhv=[28 -80.84 4.305];
% subs(8).camera_up_vector=dict();
% subs(8).camera_up_vector('topdown')=[0.08575 0.9602 0.2659];
% subs(8).camera_up_vector('back')=[0.0,-0.4,1];
% subs(8).camera_position=dict();
% subs(8).camera_position('topdown')=[-36.76 -453.5 1685];
% subs(8).camera_position('back')=[-148.64 -1633.98 -24.92];
