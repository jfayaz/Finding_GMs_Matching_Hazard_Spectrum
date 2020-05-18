## Finding GMs that match the Hazard Spectrum


Given the Time Histories, Spectra (RotD50, RotD100, etc) and dt of the Bi-Directional GMs and the Hazard Spectrum (UHS, CMS, etc) of the location, this code finds the requried number of Time Histories whose Spectrum matches closest to the Hazard Spectrum


INPUT:

Input Variables must be in the form of .mat file and must be in same directory. Input Variables of the .mat file should include:
      "acc_1"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-1 (any component of the bi-directional gm can be direction-1) {Contains GM Time-History}
      "acc_2"     -->  Cell Structure (1,n) with each cell of Matrix (1,m), where 'n' is the no. of gms and 'm' is the length of gm time-history in direction-2 (other component of the bi-directional gm) {Contains GM Time-History}
      "dt"        -->  Cell Structure (1,n) with each cell of Matrix (1,1) consisting of the dt of Time-History 
      "GM_Spectra"-->  Cell Structure (1,n) with each cell of Matrix (o,2), where 'n' is the no. of gms and 'o' is the length of spectra {Contains GM Spectra}

Other INPUTs consist of the following which are required to be entered below:
      "GM_Data_Filename" --> .mat file containing 'acc_1', 'acc_2', 'dt' and 'GM_Spectra' cell structures (as mentioned above)
      "file" --> File containing the Hazard Spectrum to match
      "sheet" --> Sheet name of the file containing the Hazard Spectrum to match
      "range" --> Excel Cell Range of the sheet containing the Hazard Spectrum to match
      "Reqd_No_of_GMs" --> Required No. of GMs for that match the spectra (if you want to match all type the total number of gms available in the 'GM_Data.mat' file)
      "Error_Type" --> '1' for 'Absolute Difference' ; '2' for 'Squared Difference; '3' for 'Cubic Difference' : For calculating error between Target and GM Spectrum 
      "Periods_to_Match_Spectra" --> Periods to Match the GM Spectra with Hazard Spectrum
      "Allowable_Scaling_Factors" --> Scaling factors for scaling Ground Motion Spectrum
      "Plot_Selected_GM_Histories" --> whether to plot Time-Histories of selected scaled Ground Motions (options: 'Yes','No')

Example Input data is given in the .mat files "GM_DATA.mat" and the "CMS.xlsx"


OUTPUT: 
Output ground motions will be provided in in two folders 'UNSCALED_GMS'and 'SCALED_GMS'. After the matching is conducted by scaling the GM spectra, the corresponding ground motion of the best spectral match will be provided in the 'SCALED_GMS' folder and the corresponding scaled version of the ground motion will be provided in the 'UNSCALED_GMS' folder with the best match scale factor.

Also the details of the match will be provided in the Excel file: 'GM_Results.xlsx'

The workspace output will be generated in the following variables:
"GM_1" & "GM_2" --> Corresponding Scaled Ground Motions of the best spectral match in directions 1 and 2
"GM_1_Unscaled" & "GM_2_Unscaled" --> Corresponding Unscaled Ground Motions of the best spectral match in directions 1 and 2
