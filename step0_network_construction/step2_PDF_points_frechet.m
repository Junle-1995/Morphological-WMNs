%% estimate the PDF with 2^4 to 2^10, calculate frechet distance, and compare
%% para
path = pwd;

Modal = {'DBM','VBM'};
Type = {'smooth','nosmooth'};
NP = 2.^(4:10);
roi_num = 48;
sub_num = 444;

Fre_dis = struct;
Fre_dis1 = struct;
results = struct;
results1 = struct;

%% PDF
for itype = 1:length(Type)
    for imodal = 1:length(Modal)
        cd(path)
        cd([Modal{imodal} '_' Type{itype}])
        data_path = pwd;

        cd(data_path)
        PDF_i = struct;

        for isub = 1:sub_num
            disp(['Now processing the data in subject (' num2str(isub) '/' num2str(sub_num)...
                ') ' Modal{imodal} '   |' datestr(clock)])

            clear Sig_roi
            load(['Signal_sub_' num2str(isub,'%04d') '.mat'])

            for iroi = 1:roi_num
                for iNP = 1:length(NP)
                    [PDF_i(iNP).PDF_Y(isub,iroi,:),PDF_i(iNP).PDF_X(isub,iroi,:)] = ...
                        gretna_PDF(Sig_roi{iroi}, NP(iNP));
                end
            end
        end

        cd(path)
        PDF_X = rmfield(PDF_i,'PDF_Y');
        save(['PDF_X_' Modal{imodal} '_' Type{itype} '.mat'],'PDF_X')
        PDF_Y = rmfield(PDF_i,'PDF_X');
        save(['PDF_Y_' Modal{imodal} '_' Type{itype} '.mat'],'PDF_Y')

        clear PDF_X PDF_Y
    end
end
clear PDF_i

%% frechet distance
cd(path)
for itype = 1:length(Type)
    for imodal = 1:length(Modal)
        load(['PDF_X_' Modal{imodal} '_' Type{itype} '.mat'])
        load(['PDF_Y_' Modal{imodal} '_' Type{itype} '.mat'])

        dis = struct;
        dis1 = struct;
        for iNP = 1:length(NP)-1
            dis(iNP).NP = [NP(iNP); NP(iNP+1)];
            dis(iNP).dis = zeros(sub_num,roi_num);
            dis1(iNP).NP = [NP(iNP); NP(iNP+1)];
            dis1(iNP).dis = zeros(sub_num,roi_num);

            for isub = 1:sub_num
                disp(['Now calculating the data in subject (' num2str(isub) '/' num2str(sub_num)...
                    ')    |' datestr(clock)])

                for iroi = 1:roi_num
                    i_PDF = [squeeze(PDF_X(iNP).PDF_X(isub,iroi,:)),squeeze(PDF_Y(iNP).PDF_Y(isub,iroi,:))];
                    j_PDF = [squeeze(PDF_X(iNP+1).PDF_X(isub,iroi,:)),squeeze(PDF_Y(iNP+1).PDF_Y(isub,iroi,:))];

                    dis(iNP).dis(isub,iroi) = frechetDistance(i_PDF, j_PDF);
                    dis1(iNP).dis(isub,iroi) = frechetDistance(i_PDF, j_PDF,[],true);
                end
            end
        end

        Fre_dis.(Modal{imodal}) = dis;
        Fre_dis1.(Modal{imodal}) = dis1;
    end
    cd(path)
    save(['PDF_fre_dis_' Type{itype} '.mat'],'Fre_dis1','Fre_dis')
end

%% compare
for itype = 1:length(Type)
    for imodal = 1:length(Modal)
        load(['PDF_fre_dis_' Type{itype} '.mat'],'Fre_dis1','Fre_dis')

        dis = zeros(roi_num,length(NP)-1);
        dis1 = zeros(roi_num,length(NP)-1);

        for iNP = 1:length(NP)-1
            dis(:,iNP) = mean(Fre_dis.(Modal{imodal})(iNP).dis,1)';
            dis1(:,iNP) = mean(Fre_dis1.(Modal{imodal})(iNP).dis,1)';
        end

        results(imodal).modal = Modal{imodal};
        results(imodal).dis = dis;
        for i = 1:size(dis,2)-1
            [T, P] = gretna_permutation_ttest_RM(dis(:,i), dis(:,i+1), per_time);
            results(imodal).T(i,1) = T.real;
            results(imodal).P(i,1) = P;
        end

        results1(imodal).modal = Modal{imodal};
        results1(imodal).dis = dis1;
        for i = 1:size(dis1,2)-1
            [T, P] = gretna_permutation_ttest_RM(dis1(:,i), dis1(:,i+1), per_time);
            results1(imodal).T(i,1) = T.real;
            results1(imodal).P(i,1) = P;
        end
    end
    cd(path)
    save(['PDF_dis_compare_' Type{itype} '.mat'],'results1','results')
end