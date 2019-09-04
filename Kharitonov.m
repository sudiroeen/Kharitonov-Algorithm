%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Copyright By:
% [1] Sudiro NIM: 16/399920/TK/44934
% at Sudiro@mail.ugm.ac.id
%
% [2] Nicholaus Valentino Dwicahyo NIM: 16/399912/TK/44926
% at nicovalentdwicahyo@gmail.com
%
% [2] Nathanael M. Tedjo Kurniawan NIM: 16/394955/TK/44247
% at nathanaelmic@gmail.com
%
% [2] Julio Marcelino Bagaskara NIM: 16/394946/TK/44238
% at julio.marcelino.bagaskara@mail.ugm.ac.id
%
% [3] Wardaniawan NIM: 16/399924/TK/44938
% at Daniawan90@gmail.com
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
clc
close all
% [an-,an+]s^n + [,]s^(n-1) +[,]s^(n-2)....+[a1-,a1+]s^1 + [a0-,a0+]s^0
%P = [60,30; 59,29; 58,28; 57,27; 56,26; 55,25; 54,24; 53,23; 52,22; 51,21; 50,20; 49,19; 48,18; 47,17; 46,16; 45,15; 44,14; 43,13; 42, 12; 41,11; 40,10; 39,9; 38,8; 37,7; 36,6; 35,5; 34,4; 33,3; 32,2; 31,1; 30,0 ];
 P = [120,60;118,58;116,56;114,54;112,52;110,50;108,48;106,46;104,44;102,42;100,40;98,38;96,36;94,34;92,32;90,30;88,28;86,26;84,24;82,22;80,20;78,18;76,16;74,14;72,12;70,10;68,8;66,6;64,4;62,2;60,0]
[baris, kolom] = size(P);
c = 0;
d = 0;
orde = (baris-1);
for row = baris: -1: 1
zz = baris - row + 1;
if mod(zz,2) == 1
p1(zz,1) = P(row,c + 1);
p2(zz,1) = P(row,~c + 1);
c = ~c;
p3(zz,1) = P(row,~d + 1);
p4(zz,1) = P(row,d + 1);
else
p1(zz,1) = P(row,~c + 1);
p2(zz,1) = P(row,c + 1);
p3(zz,1) = P(row,d + 1);
p4(zz,1) = P(row,~d + 1);
d = ~d;
end
end
if mod((orde+1),2) == 1
p1 = [p1;0]';
p2 = [p2;0]';
p3 = [p3;0]';
p4 = [p4;0]';
MatRH1 = zeros(orde+1,round(orde/2 + 1));
MatRH2 = zeros(orde+1,round(orde/2 + 1));
MatRH3 = zeros(orde+1,round(orde/2 + 1));
MatRH4 = zeros(orde+1,round(orde/2 + 1));
else
p1 = p1'
p2 = p2';
p3 = p3';
p4 = p4';
MatRH1 = zeros(orde+1,round((orde+1)/2));
MatRH2 = zeros(orde+1,round((orde+1)/2));
MatRH3 = zeros(orde+1,round((orde+1)/2));
MatRH4 = zeros(orde+1,round((orde+1)/2));
end
p1 = p1
p2 = p2
p3 = p3
p4 = p4
[barisMat, kolomMat] = size(MatRH1);
for row = barisMat:-1:1
for col = kolomMat:-1:1
m = barisMat - row + 1;
n = kolomMat - col + 1;
if m == 1
MatRH1(m,n) = p1(1,2*(kolomMat - col)+1);
MatRH2(m,n) = p2(1,2*(kolomMat - col)+1);
MatRH3(m,n) = p3(1,2*(kolomMat - col)+1);
MatRH4(m,n) = p4(1,2*(kolomMat - col)+1);
elseif m == 2
MatRH1(m,n) = p1(1,2*(kolomMat - col)+2);
MatRH2(m,n) = p2(1,2*(kolomMat - col)+2);
MatRH3(m,n) = p3(1,2*(kolomMat - col)+2);
MatRH4(m,n) = p4(1,2*(kolomMat - col)+2);
end
end
end
[barisMat, kolomMat] = size(MatRH1);
for row = barisMat:-1:2
for col = kolomMat:-1:1
m = barisMat - row + 3;
n = kolomMat - col + 1;
if n == kolomMat
MatRH1((m-2):(m-1),n+1) = zeros(2,1);
MatRH2((m-2):(m-1),n+1) = zeros(2,1);
MatRH3((m-2):(m-1),n+1) = zeros(2,1);
MatRH4((m-2):(m-1),n+1) = zeros(2,1);
end
% Antisipasi saat semua barisnya nol, sehingga diisi dengan
% auxilary equation
cekAux1 = MatRH1(m-1,1:kolomMat);
cekAux2 = MatRH2(m-1,1:kolomMat);
cekAux3 = MatRH3(m-1,1:kolomMat);
cekAux4 = MatRH4(m-1,1:kolomMat);
if cekAux1 == zeros(1,kolomMat)
pp1 = MatRH1(m-2,1:kolomMat);
for Rorde = 1:kolomMat
pp1(1,Rorde) = (row-1-2*(Rorde-1))*pp1(1,Rorde);
end
for col1 = kolomMat:-1:1
MatRH1(m-1,kolomMat-col1+1) = pp1(1,kolomMat - col1+1);
end
end
if cekAux2 == zeros(1,kolomMat)
pp2 = MatRH2(m-2,1:kolomMat);
for Rorde = 1:kolomMat
pp2(1,Rorde) = (row-1-2*(Rorde-1))*pp2(1,Rorde);
end
for col2 = kolomMat:-1:1
MatRH2(m-1,kolomMat-col2+1) = pp2(1,kolomMat - col2+1);
end
end
if cekAux3 == zeros(1,kolomMat)
pp3 = MatRH3(m-2,1:kolomMat);
for Rorde = 1:kolomMat
pp3(1,Rorde) = (row-1-2*(Rorde-1))*pp3(1,Rorde);
end
for col3 = kolomMat:-1:1
MatRH3(m-1,kolomMat-col3+1) = pp3(1,kolomMat - col3+1);
end
end
if cekAux4 == zeros(1,kolomMat)
pp4 = MatRH4(m-2,1:kolomMat);
for Rorde = 1:kolomMat
pp4(1,Rorde) = (row-1-2*(Rorde-1))*pp4(1,Rorde);
end
for col4 = kolomMat:-1:1
MatRH4(m-1,kolomMat-col4+1) = pp4(1,kolomMat - col4+1);
end
end
% Antisipasi saat kolom pertama nol, tapi tidak seluruh barisnya nol
% diberi nilai epsilon (angka positif yang sangat kecil)
% di sini diberi nilai 1e-19
format long
for r = 2:barisMat
if MatRH1(r,1) == 0
MatRH1(r,1) = 1e-19;
end
if MatRH2(r,1)== 0
MatRH2(r,1) = 1e-19;
end
if MatRH3(r,1)== 0
MatRH3(r,1) = 1e-19;
end
if MatRH4(r,1)== 0
MatRH4(r,1) = 1e-19;
end
end
if m <= barisMat
Mat1 = [MatRH1((m-2):(m-1),1),MatRH1((m-2):(m-1),n+1)];
Mat2 = [MatRH2((m-2):(m-1),1),MatRH2((m-2):(m-1),n+1)];
Mat3 = [MatRH3((m-2):(m-1),1),MatRH3((m-2):(m-1),n+1)];
Mat4 = [MatRH4((m-2):(m-1),1),MatRH4((m-2):(m-1),n+1)];
MatRH1(m,n)= - det(Mat1)/MatRH1(m-1,1);
MatRH2(m,n)= - det(Mat2)/MatRH2(m-1,1);
MatRH3(m,n)= - det(Mat3)/MatRH3(m-1,1);
MatRH4(m,n)= - det(Mat4)/MatRH4(m-1,1);
end
end
end
MatRH1 = MatRH1(1:barisMat,1:kolomMat)
MatRH2 = MatRH2(1:barisMat,1:kolomMat)
MatRH3 = MatRH3(1:barisMat,1:kolomMat)
MatRH4 = MatRH4(1:barisMat,1:kolomMat)
x = 1:barisMat;
vectorK1 = [MatRH1(x,1);MatRH2(x,1);MatRH3(x,1);MatRH4(x,1)];
[barisVek,kolomVek]=size(vectorK1);
counterTrue = 0;
for k = 1:1:barisVek
if vectorK1(k) < 0
counterTrue = counterTrue + 1;
end
end
if counterTrue > 0
counterTrue = counterTrue
disp('System is not stable')
else
disp('System is stable')
end
