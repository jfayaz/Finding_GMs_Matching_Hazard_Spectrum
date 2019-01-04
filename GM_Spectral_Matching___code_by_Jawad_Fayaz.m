clear ; close all; clc; 
%% ====================== SPECTRAL MATCHING ============================ %%
%  author : JAWAD FAYAZ (email: jfayaz@uci.edu)

%  ------------- Instructions -------------- %
%  Given the Time Histories, Spectra (RotD50, RotD100, etc) and dt of the Bi-Directional GMs 
%  and the Hazard Spectrum (UHS, CMS, etc) of the location, this code finds
%  the requried number of Time Histories whose Spectrum matches closest to 
%  the Hazard Spectrum
%  
%  INPUT:
%  Input Variables must be in the form of .mat file and must be in same directory
%  Input Variables of the .mat file should include:
%  "acc_1"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-1 (any component of the bi-directional gm can be direction-1) {Contains GM Time-History}
%  "acc_2"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-2 (other component of the bi-directional gm) {Contains GM Time-History}
%  "dt"        -->  Cell Structure (1,n) with each cell of Matrix (1,1) consisting of the dt of Time-History 
%  "GM_Spectra"-->  Cell Structure (1,n) with each cell of Matrix (o,2), where 'n' is the no. of gms and 'o' is the length of spectra {Contains GM Spectra}
% 
%  Other INPUTs consist of the following which are required to be entered below:
%  "GM_Data_Filename" --> .mat file containing 'acc_1', 'acc_2', 'dt' and 'GM_Spectra' cell structures (as mentioned above)
%  "file" --> File containing the Hazard Spectrum to match
%  "sheet" --> Sheet name of the file containing the Hazard Spectrum to match
%  "range" --> Excel Cell Range of the sheet containing the Hazard Spectrum to match
%  "Reqd_No_of_GMs" --> Required No. of GMs for that match the spectra (if you want to match all type the total number of gms available in the 'GM_Data.mat' file)
%  "Error_Type" --> '1' for 'Absolute Difference' ; '2' for 'Squared Difference; '3' for 'Cubic Difference' : For calculating error between Target and GM Spectrum 
%  "Periods_to_Match_Spectra" --> Periods to Match the GM Spectra with Hazard Spectrum
%  "Allowable_Scaling_Factors" --> Scaling factors for scaling Ground Motion Spectrum
%  "Plot_Selected_GM_Histories" --> whether to plot Time-Histories of selected scaled Ground Motions (options: 'Yes','No')
%
%  Example Input data is given in the .mat files "GM_DATA.mat" and the "CMS.xlsx"
%
%
%  OUTPUT:
%  Output ground motions will be provided in in two folders 'UNSCALED_GMS'
%  and 'SCALED_GMS' . After the matching is conducted by scaling the
%  GM spectra, the corresponding ground motion of the best spectral match
%  will be provided in the 'SCALED_GMS' folder and the corresponding scaled
%  version of the ground motion will be provided in the 'UNSCALED_GMS'
%  folder with the best match scale factor.
%  Also the details of the match will be provided in the Excel file: 'GM_Results.xlsx'
%  
%  The workspace output will be generated in the following variables:
%  "GM_1" & "GM_2" --> Corresponding Scaled Ground Motions of the best spectral match in directions 1 and 2
%  "GM_1_Unscaled" & "GM_2_Unscaled" --> Corresponding Unscaled Ground Motions of the best spectral match in directions 1 and 2


%% ========================= USER INPUTS =============================== %%

GM_Data_Filename = 'GM_DATA.mat';               % .mat file containing 'acc_1', 'acc_2', 'dt' and 'GM_Spectra' cell structures (as mentioned above)
file             = 'CMS.xlsx';                  % File containing the Hazard Spectrum to match
sheet            = 'Sheet1';                    % Sheet name of the file containing the Hazard Spectrum to match
range            = 'A2:E100';                   % Excel Cell Range of the sheet containing the Hazard Spectrum to match
Reqd_No_of_GMs   = 2;                           % Required No. of GMs for that match the spectra (if you want to match all type the total number of gms available in the 'GM_Data.mat' file)
Error_Type       = 2;                           % '1' for 'Absolute Difference' ; '2' for 'Squared Difference; '3' for 'Cubic Difference' : For calculating error between Target and GM Spectrum 
Periods_to_Match_Spectra  = [0.1:0.1:2.0];      % Periods to Match the GM Spectra with Hazard Spectrum
Allowable_Scaling_Factors = [0.2:0.2:10];       % Scaling factors for scaling Ground Motion Spectrum
Plot_Selected_GM_Histories= 'Yes';              % Plot Time-Histories of selected scaled Ground Motions (options: 'Yes','No')

%%%%%%================= END OF USER INPUT ========================%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------- CALLING FUNCTIONS ----------------

load(GM_Data_Filename)
if Reqd_No_of_GMs > length(GM_Spectra)
    fprintf('Reqd_No_of_GMs cannot be greater than number of ground motions provided in the "GM_Data.mat" file')
end

[Interpolated_Spectrum,User_Defined_Spectrum,GM_Spectra] = interpole_ars_spectrum(file,sheet,range,GM_Spectra,Periods_to_Match_Spectra);

fprintf('Scaling and Checking Differences between Target and Ground Motion Spectrum...\n')
[Stepped_Error_Rec,Error_Matrix] = diff_bw_GivenSpectrum_GMSpectrum(Periods_to_Match_Spectra,GM_Spectra,Interpolated_Spectrum,Allowable_Scaling_Factors,Error_Type);
[Error_Matrix_for_Each_GM,Min_Error_for_Each_GM] = arrange_errors(GM_Spectra,Error_Matrix,Allowable_Scaling_Factors);
    
fprintf('Finding Match for the Target Spectra.....\n')
[Selected_GMs_Data,Closest_Scaled_Spectra,Closest_Unscaled_Spectra,Sorted_Min_Error_for_Each_GM] = find_matching_spectras(Interpolated_Spectrum,User_Defined_Spectrum,Min_Error_for_Each_GM,Reqd_No_of_GMs,GM_Spectra);
[GM_1_Unscaled, GM_1, GM_2_Unscaled, GM_2, Delta_T, NPTS] = derive_GMs_matching_spectra(Reqd_No_of_GMs,Selected_GMs_Data,acc_1,acc_2,dt,Plot_Selected_GM_Histories);

fprintf('Writing Ground Motions to .AT2 files in the current directory.....\n')
write_GMs_to_AT2(Selected_GMs_Data,Reqd_No_of_GMs,GM_1,GM_2,Delta_T,GM_1_Unscaled,GM_2_Unscaled);

fprintf('Writing GM Spectra to XLSX File.....\n')
[Selected_GMs_Data, Sorted_Min_Error_for_Each_GM] = write_spectra_to_excel(GM_Spectra,Reqd_No_of_GMs,Selected_GMs_Data,NPTS,Delta_T,Closest_Scaled_Spectra,Closest_Unscaled_Spectra,Sorted_Min_Error_for_Each_GM);
    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% --------------- INTERPOLATION OF ARS SPECTRUM ---------------
function [Interpolated_Spectrum,User_Defined_Spectrum,GM_Spectra] = interpole_ars_spectrum(file,sheet,range,GM_Spectra,Periods_to_Match_Spectra)
    T = Periods_to_Match_Spectra ;
    ars_data = xlsread(file,sheet,range);
    User_Defined_Spectrum = [ars_data(:,1),ars_data(:,2)];
    Interpolated_Sa = interp1(User_Defined_Spectrum(:,1),User_Defined_Spectrum(:,2),T');
    Interpolated_Spectrum = [T' ,  Interpolated_Sa];
    [r,~] = find(isnan(Interpolated_Spectrum)==1);
    Interpolated_Spectrum(r,:) = [];

    scatter (User_Defined_Spectrum(:,1),User_Defined_Spectrum(:,2),'r','LineWidth',2);
    title ('Hazard Spectrum','fontsize',16,'fontWeight','bold')
    xlabel ('Period (sec)','fontsize',16,'fontWeight','bold')
    ylabel ('Spectral Acceleration (g)','fontsize',16,'fontWeight','bold')
    hold on
    plot(Interpolated_Spectrum(:,1),Interpolated_Spectrum(:,2),'b*-','LineWidth',2)
    set(gca,'fontsize',14,'FontName', 'Times New Roman','LineWidth', 1.5)
    grid on; box off
    legend ('User-Defined Spectrum','Interpolated Portion for Matching')
    
    for i = 1:length(GM_Spectra)
        GM_Spectra{i}(r,:) = [];
    end
end


%% ----- Difference between Interpolated ARS Spectrum and Ground Motion Spectrum ----- %%
function [Stepped_Error_Rec,Error_Matrix] = diff_bw_GivenSpectrum_GMSpectrum(Periods_to_Match_Spectra,GM_Spectra,Interpolated_Spectrum,Allowable_Scaling_Factors,Error_Type)
    T = Periods_to_Match_Spectra;
    Area_ARS = trapz(Interpolated_Spectrum(:,1), Interpolated_Spectrum(:,2));
    Area_Cum_ARS = cumtrapz(Interpolated_Spectrum(:,1), Interpolated_Spectrum(:,2));

    for v = 1:(length(Area_Cum_ARS)-1)
        Area_Stepped_ARS(v) = Area_Cum_ARS(v+1)-Area_Cum_ARS(v);
    end

    d = 0;
    for sf = 1:length(Allowable_Scaling_Factors)
        c = 0;
        fprintf('Scaling and Matching GM Spectra scaled with %.2f\n',Allowable_Scaling_Factors(sf))
        for gm = 1:length(GM_Spectra)       
            c = c +1;
            d = d+1;
            data = GM_Spectra{gm};
            data_sf(:,1) = data(:,1);
            data_sf(:,2) = data(:,2).*Allowable_Scaling_Factors(sf); 
            Matching_nos = interp1(data(:,1), data(:,1), T, 'nearest');
            [~,idx] = ismember(data(:,1),Matching_nos);
            rows =  find(idx > 0);
            Area_GM = trapz(data_sf(rows,1), data_sf(rows,2));
            Area_Cum_GM = cumtrapz (data_sf(rows,1), data_sf(rows,2));
            for w = 1:(length(Area_Cum_GM)-1)
                Area__Stepped_GM(w) = Area_Cum_GM(w+1)-Area_Cum_GM(w);
            end
            Stepped_Error = (Area__Stepped_GM - Area_Stepped_ARS).^Error_Type;
            Stepped_Error_Rec{sf}{gm} = abs(Stepped_Error);
            Error_Matrix{sf}{gm} =  sum(abs(Stepped_Error));
            Error_Array(d,1) = sum(abs(Stepped_Error));
            clearvars data 
        end   
    end
end
    

%% --------------- Arranging Errors as per Realisations --------------- %%
function [Error_Matrix_for_Each_GM,Min_Error_for_Each_GM] = arrange_errors(GM_Spectra,Error_Matrix,Allowable_Scaling_Factors)
    for i = 1:length(GM_Spectra)
        for j = 1:length(Allowable_Scaling_Factors)
            Error_Matrix_for_Each_GM {i}(:,j) = Error_Matrix{j}{i};
        end
    end

    k = 0;
    for i = 1:length(Error_Matrix_for_Each_GM)
        k = k+1;
        [~, n] = find(Error_Matrix_for_Each_GM{i} == min(Error_Matrix_for_Each_GM{i}));
        Min_Error_for_Each_GM(k,1) = min(Error_Matrix_for_Each_GM{i});                %% Minimum Error for each ground motion among all scaling factors
        Min_Error_for_Each_GM(k,2) = i;                                               %% Ground Motion Number                                                            
        Min_Error_for_Each_GM(k,3) = Allowable_Scaling_Factors(n);                    %% Scaling Factor corresponding to Minimum Error of the corresponding Ground Motion
    end
end


%% --------------- Matching Spectras --------------- %%
function [Selected_GMs_Data,Closest_Scaled_Spectra,Closest_Unscaled_Spectra,Sorted_Min_Error_for_Each_GM] = find_matching_spectras(Interpolated_Spectrum,User_Defined_Spectrum,Min_Error_for_Each_GM,Reqd_No_of_GMs,GM_Spectra)
    [~, index] = sort(Min_Error_for_Each_GM(:,1));
    Sorted_Min_Error_for_Each_GM = Min_Error_for_Each_GM(index,:);
    
    figure()
    scatter (User_Defined_Spectrum(:,1),User_Defined_Spectrum(:,2),80,'b','filled');
    xlabel ('Period (sec)')
    ylabel ('Spectral Acceleration (g)')
    hold on
    plot(Interpolated_Spectrum(:,1),Interpolated_Spectrum(:,2),'r','LineWidth',3);
    legend('ARS Spectrum','Interpolated Spectrum') 
    
    for i =1:Reqd_No_of_GMs
        GM_Num = Sorted_Min_Error_for_Each_GM (i,2);
        GM_Sf = Sorted_Min_Error_for_Each_GM (i,3);
        Selected_GMs_Data(i,1) = GM_Num;
        Selected_GMs_Data(i,2) = GM_Sf;
        spectr  =  GM_Spectra{GM_Num};
        Sc_Spectrum(:,1) = spectr (:,1);
        Sc_Spectrum(:,2) = spectr (:,2).*GM_Sf;
        Un_Sc_Spectrum(:,1) = spectr (:,1);
        Un_Sc_Spectrum(:,2) = spectr (:,2);
        Closest_Scaled_Spectra{i} = cell2struct((num2cell(Sc_Spectrum))',{'Period','RotD50'});
        Closest_Unscaled_Spectra{i} = cell2struct((num2cell(Un_Sc_Spectrum))',{'Period','RotD50'});
        
        if i > 1
            plot(Sc_Spectrum(:,1),Sc_Spectrum(:,2),'-','Color',[0.4 0.4 0.4],'LineWidth',1.5,'HandleVisibility','off');          
        else
            plot(Sc_Spectrum(:,1),Sc_Spectrum(:,2),'-','Color',[0.4 0.4 0.4],'LineWidth',1.5,'DisplayName',[num2str(Reqd_No_of_GMs),' Closest Matching GM Spectra']);          
        end
        
        clearvars spectrum
    end
end


%% --------------- Deriving Ground Motions Matching the Spectra --------------- %%
function [GM_1_Unscaled, GM_1, GM_2_Unscaled, GM_2, Delta_T, NPTS] = derive_GMs_matching_spectra(Reqd_No_of_GMs,Selected_GMs_Data,acc_1,acc_2,dt,Plot_Selected_GM_Histories)
    figure('Name','SELECTED GROUND MOTIONS` TIME HISTORIES')
    for j = 1:size(Selected_GMs_Data,1)
        GM_1_Unscaled{j} = acc_1{Selected_GMs_Data(j,1)};
        GM_1{j} = GM_1_Unscaled{j}.*Selected_GMs_Data(j,2);
        GM_2_Unscaled{j} = acc_2{Selected_GMs_Data(j,1)};
        GM_2{j} = GM_2_Unscaled{j}.*Selected_GMs_Data(j,2);
        Delta_T(j) = dt{Selected_GMs_Data(j,1)};
        NPTS(j) = length(GM_1{j});
        
        if strcmpi(Plot_Selected_GM_Histories,'Yes') == 1
            subplot(ceil(sqrt(Reqd_No_of_GMs/2)),2*ceil(sqrt(Reqd_No_of_GMs/2)),j)
            if j > 1
                plot([Delta_T(1):Delta_T(1):(Delta_T(1)*length(GM_1{j}))],GM_1{j},'r','Linewidth',1.25,'HandleVisibility','off')
                hold on
                plot([Delta_T(1):Delta_T(1):(Delta_T(1)*length(GM_2{j}))],GM_2{j},'b','Linewidth',1.25,'HandleVisibility','off')
            else
                plot([Delta_T(1):Delta_T(1):(Delta_T(1)*length(GM_1{j}))],GM_1{j},'r','Linewidth',1.25,'DisplayName','Time History in Dir 1')
                hold on
                plot([Delta_T(1):Delta_T(1):(Delta_T(1)*length(GM_2{j}))],GM_2{j},'b','Linewidth',1.25,'DisplayName','Time History in Dir 2')
                lgd = legend();
                lgd.FontSize = 14;
            end
            xlabel ('Time(s)','fontWeight','bold')
            ylabel ('Acc(g)','fontWeight','bold')
        end
    end
end


%% --------------- Writing Ground Motions to .AT2 file --------------- %%
function write_GMs_to_AT2(Selected_GMs_Data,Reqd_No_of_GMs,GM_1,GM_2,Delta_T,GM_1_Unscaled,GM_2_Unscaled)
    directory = pwd;
    mkdir SCALED_GMS
    cd ([directory,'\SCALED_GMS']);
    for l = 1:Reqd_No_of_GMs
        %%%% Writing GM history in direction 1 (North) %%%%
        Direction = ['GM1',num2str(l),'.AT2'];
        gm1 = GM_1{l}';
        npts = length(gm1);
        fid = fopen (Direction,'w');
        fprintf(fid,'Closest Spectral Match Scaled Ground Motion in Direction-1 (units g)\n');
        fprintf(fid,'This Ground Motion is the Ground Motion No. %.0f of the User-Provided Database\n',Selected_GMs_Data(l,1));
        fprintf(fid,'This Ground Motion is the scaled version with a Scale Factor = %.3f\n',Selected_GMs_Data(l,2));
        fprintf(fid,'NPTS=  %d, DT= %f\n',npts,Delta_T(l));
        fprintf (fid,'%d\n',gm1);
        fclose(fid);  

        %%%% Writing GM history in direction 2 (East) %%%%
        Direction = ['GM2',num2str(l),'.AT2'];
        gm2 = GM_2{l}';
        npts = length(gm2);
        fid = fopen (Direction,'w');
        fprintf(fid,'Closest Spectral Match Scaled Ground Motion in Direction-2 (units g)\n');
        fprintf(fid,'This Ground Motion is the Ground Motion No. %.0f of the User-Provided Database\n',Selected_GMs_Data(l,1));
        fprintf(fid,'This Ground Motion is the scaled version with a Scale Factor = %.3f\n',Selected_GMs_Data(l,2));        fprintf(fid,'NPTS=  %d, DT= %f\n',npts,Delta_T(l));
        fprintf (fid,'%d\n',gm2);
        fclose(fid); 

        clearvars gm1 gm2
    end

    cd (directory)

    mkdir UNSCALED_GMS
    cd ([directory,'\UNSCALED_GMS'])
    for l = 1:Reqd_No_of_GMs
        %%%% Writing GM history in direction 1 (North) %%%%
        Direction = ['GM1',num2str(l),'.AT2'];
        gm1 = GM_1_Unscaled{l}';
        npts = length(gm1);
        fid = fopen (Direction,'w');
        fprintf(fid,'Closest Spectral Match corresponding Unscaled Ground Motion in Direction-1 (units g)\n');
        fprintf(fid,'This Ground Motion is the Ground Motion No. %.0f of the User-Provided Database\n',Selected_GMs_Data(l,1));
        fprintf(fid,'This Ground Motion is the unscaled version hence with a Scale Factor = 1\n');
        fprintf(fid,'NPTS=  %d, DT= %f\n',npts,Delta_T(l));
        fprintf (fid,'%d\n',gm1);
        fclose(fid);  

        %%%% Writing GM history in direction 2 (East) %%%%
        Direction = ['GM2',num2str(l),'.AT2'];
        gm2 = GM_2_Unscaled{l}';
        npts = length(gm2);
        fid = fopen (Direction,'w');
        fprintf(fid,'Closest Spectral Match corresponding Unscaled Ground Motion in Direction-2 (units g)\n');
        fprintf(fid,'This Ground Motion is the Ground Motion No. %.0f of the User-Provided Database\n',Selected_GMs_Data(l,1));
        fprintf(fid,'This Ground Motion is the unscaled version hence with a Scale Factor = 1\n');
        fprintf(fid,'NPTS=  %d, DT= %f\n',npts,Delta_T(l));
        fprintf (fid,'%d\n',gm2);
        fclose(fid); 

        clearvars gm1 gm2
    end
    cd (directory)
end


%% --------------- Generating XLSX file --------------- %%
function [Selected_GMs_Data, Sorted_Min_Error_for_Each_GM] = write_spectra_to_excel(GM_Spectra,Reqd_No_of_GMs,Selected_GMs_Data,NPTS,Delta_T,Closest_Scaled_Spectra,Closest_Unscaled_Spectra,Sorted_Min_Error_for_Each_GM)  
    writefile = 'GM_Results.xlsx';
    headings = strsplit([ 'GM_Name_X ', 'GM_Name_Y ', 'Ground_Motion_No ' , 'NPTS ', 'dt ' , 'Scale_Factor'],' ');
    title1{1} = ['GROUND_MOTION_SCALED_SPECTRA_RotD50 (g)'];
    title2{1} = ['GROUND_MOTION_UNSCALED_SPECTRA_RotD50 (g)'];
    GM_Names_Array{1,1} = 'Tn (sec)';

    for m = 1:Reqd_No_of_GMs
        GM_Names{m,1} = ['GM1',num2str(m)];
        GM_Names{m,2} = ['GM2',num2str(m)];
        GM_Names_Array{m+1,1} = ['GM_',num2str(m)];
    end

    write = [Selected_GMs_Data(:,1) , NPTS' , Delta_T' ,Selected_GMs_Data(:,2)]; 

    xlswrite (writefile,headings,'A1:F1')
    xlswrite (writefile,write, ['C2:F',num2str(2+size(write,1)-1)])
    xlswrite (writefile,GM_Names,['A2:B',num2str(2+size(write,1)-1)])
    xlswrite (writefile,title1,['A',num2str(7+size(write,1)-1),':A',num2str(7+size(write,1)-1)])
    xlswrite (writefile,title2,['A',num2str(12 + 2*size(write,1)-1),':A',num2str(12 + 2*size(write,1)-1)])

    TN = GM_Spectra{1}(:,1)';
    SCALE_spectrum(:,1) = TN';
    UNSCALE_spectrum(:,1) = TN';
    for n = 1:length(Closest_Scaled_Spectra)
        sc_spctrm  =  cell2mat(struct2cell(Closest_Scaled_Spectra{n})');
        un_sc_spctrm  =  cell2mat(struct2cell(Closest_Unscaled_Spectra{n})');
        SCALE_spectrum(:,n+1) = sc_spctrm (:,2);
        UNSCALE_spectrum(:,n+1) = un_sc_spctrm (:,2);
        clearvars un_sc_spctrm sc_spctrm
    end

    %%% Writing for Scaled Ground Motions
    xlswrite (writefile,GM_Names_Array,['A',num2str(7+size(write,1)),':A',num2str(7+2*size(write,1))])
    xlswrite (writefile,SCALE_spectrum',['B',num2str(7+size(write,1)),':BU',num2str(7+2*size(write,1))])

    %%% Writing for UnScaled Ground Motions
    xlswrite (writefile,GM_Names_Array,['A',num2str(12+2*size(write,1)),':A',num2str(12+3*size(write,1))])
    xlswrite (writefile,UNSCALE_spectrum',['B',num2str(12+2*size(write,1)),':BU',num2str(12+3*size(write,1))])

    Selected_GMs_Data = cell2struct(num2cell(Selected_GMs_Data),{'GM_No','Scale_Factor'},2);
    Sorted_Min_Error_for_Each_GM = cell2struct(num2cell(Sorted_Min_Error_for_Each_GM),{'Error','GM_No','Scale_Factor'},2);

    clearvars -except acc_2 acc_1 simulations_scenarios is_pulse dt Closest_Scaled_Spectra Delta_T Error_Matrix GM_1 GM_2 Interpolated_ARS_Spectrum NPTS Reqd_No_of_GMs Allowable_Scaling_Factors GM_Spectra Stepped_Error_Rec TN file sheet range Closest_Unscaled_Spectra GM_2_Unscaled GM_1_Unscaled Sorted_Min_Error_for_Each_GM Selected_GMs_Data T
end
                    
 