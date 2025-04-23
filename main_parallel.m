%villyes, 20.04.25

% Example call for hex_class-call; for PCHE-type heat exchangers
% make sure to install Parallel Computing Toolbox
%% INPUTS: example data for recalculation, SI-units
% if using experimental data and later comparing, make sure to make the
% data fit the energy conversion on both sides exactly

tab=table;
tab.id=["dataset_1"; "dataset_2"];
tab.mCO2=[0.6;0.4]; %kg/s (only given once, since recuperator.)
tab.THotIn=[360;450]+273.15; %K
tab.pHotIn=[61*10^5;60*10^5]; %Pa (absolute pressure, if using experimental data make sure to add ambient p)
tab.TColdIn=[20;20]+273.15; %K
tab.pColdIn=[150*10^5;200*10^5]; %Pa

%% INPUTS: prepare objects with geometrical data
pche_zigzag;
HTXzigzag.nCells=200; %discretization

n_r=height(tab);

for i=1:n_r
    HTXa(i,1)=copy(HTXzigzag);
end

clear HTXzigzag
%% give data to recalculate

for i=1:n_r
    HTXa(i,1).name = tab.id(i);
    HTXa(i,1).Qdot_start = 2*10^5; %if you have values, use them here
    HTXa(i,1).T1_start = tab.THotIn(i);
    HTXa(i,1).T3_start = tab.TColdIn(i);
    HTXa(i,1).mDot1 = tab.mCO2(i); 
    HTXa(i,1).mDot3 = tab.mCO2(i);
    HTXa(i,1).p1in = tab.pHotIn(i); 
    HTXa(i,1).p3in = tab.pColdIn(i);
    HTXa(i,1).delta_p1 = 1*10^5;
    HTXa(i,1).delta_p3 = 0.5*10^5;
    HTXa(i,1).HTcorr = "dittusboelter";
end


%% feed correlations
HTcorr=["meshram","meshram2","dittusboelter","gnielinksi_konakov","kim","saeedkim","cheng"]; %"overallHTC", "meshram", "meshram2", "kim", "saeedkim", "cheng", "dittusboelter", "jacksonpitla", "gnielinski_filonenko", "gnielinski_konakov"
n_c=length(HTcorr);
HTXa=repmat(HTXa,1,n_c);

for i=1:n_c
    for j=1:n_r
        HTXa(j,i)=copy(HTXa(j,1));
        HTXa(j,i).HTcorr=HTcorr(i);
    end        
end

%%  calulating! (1st column) 
% calculations are done column-wise (all datasets with one heat
% correlation). Could also be done element-wise but saving only occurs at
% the end.

for i=1:n_r %non-parallel for i=1:n_r   %element-wise parfor i=1:n_r*n_c
    local=HTXa(i,1);  %elementwise local=HTXa(i);
    tic
    %if local.calculated==0
        TQ_diagram(local)

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,1)=local;
    %end
end
disp('Column 1 calculated.')
save('save_hexclass_pche')

%  calulating! (2nd column) 
parfor i=1:n_r
    local=HTXa(i,1);
    tic
    %if local.calculated==0
        TQ_diagram(local)

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,1)=local;
    %end
end
disp('Column 2 calculated.')
% save('save_hexclass_pche')

% calulating! (3rd column)
parfor i=1:n_r
    local=HTXa(i,3);
    tic
    %if local.calculated==0
        TQ_diagram(local);

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,3)=local;
    %end
end
disp('Column 3 calculated.')
% save('save_hexclass_pche')

% calulating! (4th column)
parfor i=1:n_r
    local=HTXa(i,4);
    tic
    %if local.calculated==0
        TQ_diagram(local);

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,4)=local;
    %end
end
disp('Column 4 calculated')
% save('save_hexclass_pche')

% calulating! (5th column)
parfor i=1:n_r
    local=HTXa(i,5);
    tic
    %if local.calculated==0
        TQ_diagram(local);

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,5)=local;
    %end
end
disp('Column 5 calculated.')
% save('save_hexclass_pche')

% calulating! (6th column)
parfor i=1:n_r
    local=HTXa(i,6);
    tic
    %if local.calculated==0
        TQ_diagram(local);

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,6)=local;
    %end
end
disp('Column 6 calculated.')
% save('save_hexclass_pche1')

% calulating! (7th column)
parfor i=1:n_r
    local=HTXa(i,7);
    tic
    %if local.calculated==0
        TQ_diagram(local);

        local.hex_calculation(local);
        local.calctime=toc;

        HTXa(i,7)=local;
    %end
end
disp('Column 7 calculated.')
% save('save_hexclass_pche1')
