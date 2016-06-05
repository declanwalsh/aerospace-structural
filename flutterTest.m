% short script to test the flutter calc function

close all
clear all

c = 2.5*25.4;
freqTwist = 150;
freqBend = 40;
a = -0.15;

[vFlut, vDiv] = flutterCalc(c, freqTwist, freqBend, a)