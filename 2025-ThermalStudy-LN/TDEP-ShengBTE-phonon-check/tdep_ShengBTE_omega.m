%对比ShengBTE的频率和Phonopy
clc;clear all
tdep=load('.\outfile.dispersion_relations');
m=size(tdep,1);
index=[1:2:m,m];
index2=index(2:end);%去除gamma
% 这是因为ShengBTE未对gamma点进行校正

BTE=load('.\BTE.omega');
BTE_2=BTE/2/pi;

phonon_branch=size(BTE,2);


for i=1:phonon_branch
    plot(tdep(index2',i+1),BTE_2(2:end,i),'-v')
    hold on
end

xlabel('TDEP')
ylabel('ShengBTE')