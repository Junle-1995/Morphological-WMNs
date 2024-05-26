%% calculate the heritability of morphological similarity
path = pwd;
Output_path_data = [path '\VBM\handle'];
Output_path_h2 = [path '\VBM\result'];
mkdir(Output_path_data)
mkdir(Output_path_h2)
mkdir([path '\VBN\temporary'])

%% data reorganization
load twin_edge_data_VBM.mat
Netpara =twin_data;
tmp = zeros(3,size(Netpara,1)/2,size(Netpara,2));

for i=1:217
    m=2*i-1;
    n=2*i;
    tmp(1,i,:)=Netpara(m,:);
    tmp(2,i,:)=Netpara(n,:);
end

TwinLabel = [ones(138,1);2*ones(79,1)];
tmp(3,:,:) = repmat(TwinLabel,1,size(Netpara,2));

RealData = permute(tmp,[2 1 3]);


save([Output_path_data '\tmp_data.mat'], 'RealData');

%% calculate real h2
Rpath = 'D:\Program Files\R\R-4.2.0\bin'; % D:\Program Files\R\R-4.2.0\bin 
RscriptFileName = [path '\script.R'];
RunRcode(RscriptFileName, Rpath);

Niter = 10000;
heritability = zeros(Niter+1,size(Netpara,2));

cd(Output_path_h2)
load('tmp_data.mat')
heritability(1,:) = h2;

%% calculate rand h2
for iper = 1:Niter
    disp(['Now calculating the heritability of rand (' num2str(iper) '\'...
        num2str(Niter) '  |' datestr(clock)])
    RandIndex = randperm(size(Netpara,1)/2);
    tmp(3,:,:) = repmat(TwinLabel(RandIndex),1,size(Netpara,2));
    
    RandData = permute(tmp,[2 1 3]);
    
    save([Output_path_data '\tmp_data.mat'], 'RandData');
    RunRcode(RscriptFileName, Rpath);
    
    cd(Output_path_h2)
    load('tmp_data.mat')
    heritability(iper+1,:) = h2;
    
    save('heritability_VBM.mat','heritability');
end

