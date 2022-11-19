% Function to create time series for file input
clc; 
clear all; 
close all;

omega = 1.26;    % rad/sec
te    = 200;     % t - end
t1    = 120;     % t - e-funktion
dt    = 1;       % delta time

% first time period
t_1 = [0:dt:t1];
b   = log(omega+1)/t1;
y1  = exp(b*t_1)-1;

% second time period
t_2 = [t1+dt:dt:te];
y2  = ones(length(t_2),1)*omega

figure(); hold on; grid on;
xlabel('Zeit in sec')
T = [t_1,t_2];
Y = [y1,y2'];
plot(T,Y);
ylabel('Winkelgeschwindigkeit in rad/s')

% write into text file
fid = fopen('amplitude1.txt','w');
fprintf(fid,'!! \n');
fprintf(fid,'!! \n');
fprintf(fid,'!! number of samples \n');
fprintf(fid,' %i \n',length(T));
fprintf(fid,'!! \n');
fprintf(fid,'!! \n');
fprintf(fid,'!! time (1)    amplitude (2) \n');
for i = 1:length(T)
    fprintf(fid,'%12.5f %12.5f \n',[T(i),Y(i)]);
end
fclose(fid);
return