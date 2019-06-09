clc;
clear;

load('R.mat')
load('Base1_B.mat')
snr_start = -1;
p = 2;
snr_out = RCA_threshold_predict(Base_1,snr_start,p,R,snr_R);