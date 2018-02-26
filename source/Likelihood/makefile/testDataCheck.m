%===================================
%
%   Test Synthetic data
%
%==================================



function testDataCheck

close all;
clc


addpath('/Users/lipkova/WORK/UQ/UQ_VillaGarbald/Matlab2C/matrixMatlab2Cpp/matlab/')


% A=loadMatrix('/cluster/home/mavt/lipkovaj/GliomaAdvance/source/Anatmoy/Patient01/FET1_1b.dat');
A=loadMatrix('HGG_data.dat');
size(A)
S=loadMatrix('Brain.dat');
size(S)
%%
[Nx,Ny,Nz] = size(A)
B=zeros(Nx,Ny,Nz);
S = S';

[Px,Py]=size(S);

for i = 1:Px
    ix = S(i,1);
    iy = S(i,2);
    iz = S(i,3);
    
    B(ix,iy,iz) = 10;
    
end;

fid = 1;

for iz=64:64
    figure(fid)
    pcolor(A(:,:,iz))
    fid = fid+1;
    
    figure(fid)
    pcolor(B(:,:,iz))
    fid = fid+1;
    colorbar;
end;



% figure(fid)
% pcolor(A(:,:,63))
% fid = fid+1;
%
% figure(fid)
% pcolor(A(:,:,62))
% fid = fid+1;
%
% figure(fid)
% pcolor(A(:,:,61))
% fid = fid+1;



