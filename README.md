## Finding GMs that match the Hazard Spectrum


Given the Time Histories, Spectra (RotD50, RotD100, etc) and dt of Bi-Directional GMs or using the NGAWest2 Database, this code utilizes the user-provided target Hazard Spectrum (UHS, CMS, etc) to find the requried number of Time Histories whose spectra matches closest to the target Hazard Spectrum


INPUTS:


      "file" --> File containing the Hazard Spectrum to match (example is provided in "CMS.xlsx")

      "sheet" --> Sheet name of the file containing the Hazard Spectrum to match (example is provided in "CMS.xlsx")

      "Reqd_No_of_GMs" --> Required No. of GMs for that match the spectra (if you want to match all type the total number of gms available in the 'GM_Data.mat' or 'NGASpec.mat' file)

      "Error_Type" --> '1' for 'Absolute Difference' ; '2' for 'Squared Difference; '3' for 'Cubic Difference' : For calculating error between Target and GM Spectrum 

      "Periods_to_Match_Spectra" --> Periods to Match the GM Spectra with Hazard Spectrum

      "Allowable_Scaling_Factors" --> Scaling factors for scaling Ground Motion Spectrum

      "NGAData_or_UserData" --> This specifies whether the user wants to specify their own data with accelerations or utilize NGAWest2 database (options: 'User','NGA')

            If the user specifies 'User' as the option for the "NGAData_or_UserData" variable, then they need to provide a .mat file in the same directory
            
            The name of the .mat file must be gven to the variable "GM_Data_Filename" and the file must contain:

                  "acc_1"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-1 (any component of the bi-directional gm can be direction-1) {Contains GM Time-History}
                  
                  "acc_2"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-2 (other component of the bi-directional gm) {Contains GM Time-History}

                  "dt"        -->  Cell Structure (1,n) with each cell of Matrix (1,1) consisting of the dt of Time-History 
                  
                  "GM_Spectra"-->  Cell Structure (1,n) with each cell of Matrix (o,2), where 'n' is the no. of gms and 'o' is the length of spectra {Contains GM Spectra}
            Example Input data is given in the .mat files "GM_DATA.mat". With this option the user also needs to specify their preference for:
                  "Plot_Selected_GM_Histories" --> whether to plot Time-Histories of selected scaled Ground Motions (options: 'Yes','No')
   
      If the user specifies 'NGA' as the option for the "NGAData_or_UserData" variable, then the code will utilize the provided GMSpec file from the same directory


OUTPUT: 

If the user specifies 'User' as the option for the "NGAData_or_UserData" variable, the the following outputs will be provided: 
        
        Output ground motions will be provided in in two folders 'UNSCALED_GMS' and 'SCALED_GMS' . After the matching is conducted by scaling the GM spectra, the corresponding ground motion of the best spectral match will be provided in the 'SCALED_GMS' folder and the corresponding scaled version of the ground motion will be provided in the 'UNSCALED_GMS' folder with the best match scale factor. Also the metadata of the match will be provided in the Excel file: 'GM_Results.xlsx'
        
        The workspace output will be generated in the following variables:

            "GM_1" & "GM_2" --> Corresponding Scaled Ground Motions of the best spectral match in directions 1 and 2

            "GM_1_Unscaled" & "GM_2_Unscaled" --> Corresponding Unscaled Ground Motions of the best spectral match in directions 1 and 2
       
If the user specifies 'NGA' as the option for the "NGAData_or_UserData" variable, the the following outputs will be provided:
      Metadata including RSN, EQID, Mw, Rrup, and Vs30 will be provided in the Excel file: 'GM_Results.xlsx'
      The Excel file: 'GM_Results.xlsx' will also contain the Scaled and Unscaled spectra of the selcted ground motions
      The user can the use the RSNs to download the GMs from PEER Database

