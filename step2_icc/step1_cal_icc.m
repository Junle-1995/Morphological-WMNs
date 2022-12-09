%% calculate the ICC of each connectivity in SWU data
%% ICC script downloaded from: https://ww2.mathworks.cn/matlabcentral/mlc-downloads/downloads/submissions/22122/versions/12/download/zip/IPN_agreement.zip
%% para
sub_num = 120;
ob_k = 3;
node_num = 48;

pro_type = {'No_smooth','Smooth'};
net_type = {'VBM','DBM'};

path = pwd;

results = struct;
results(1) = [];


%% cal
for inet = 1:length(net_type)
    for ipro = 1:length(pro_type)
        results(end+1).net_type = net_type{inet};
        results(end).pro_type = pro_type{ipro};
        
        % data collect
        cd(path)
        cd ..
        data_path = [pwd '\step0_network_construction\' net_type{inet} '_SWU'];
        data = zeros(node_num,node_num,sub_num*ob_k);
        for isub = 1:ob_k*sub_num
            clear SIM
            load([data_path '\' pro_type{ipro} '\JS_KSDENSITY_256_Signal_sub_' num2str(isub,'%04d') '.mat'])
            data(:,:,isub) = SIM;
        end
        
        % cal icc
        cd(path)
        cd ..
        results(end).icc_short = zeros(node_num);
        results(end).icc_long = zeros(node_num);
        results(end).icc_all = zeros(node_num);
        
        for inode = 1:node_num
            for jnode = inode+1:node_num
                data1 = reshape(squeeze(data(inode,jnode,:)),sub_num,ob_k);
                
                results(end).icc_short(inode,jnode) = IPN_icc(data1(:,[1 2]),1,'single');
                results(end).icc_long(inode,jnode) = IPN_icc(data1(:,[1 3]),1,'single');
                results(end).icc_all(inode,jnode) = IPN_icc(data1,1,'single');
            end
        end
        results(end).icc_short = results(end).icc_short + results(end).icc_short';
        results(end).icc_long = results(end).icc_long + results(end).icc_long';
        results(end).icc_all = results(end).icc_all + results(end).icc_all';
    end
end

%% save
cd(path)
save('ICC_connectivity.mat','results')