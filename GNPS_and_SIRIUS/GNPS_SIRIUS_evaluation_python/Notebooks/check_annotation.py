# Author Louis Felix Nothias Feb 2021
# Check the SIRIUS annotation against the GNPS spectral library match consistency (support either regular/analogue modes)

import numpy as np
import pandas as pd

def check_matching_annotations(table, zodiac_score_thresh=0.7, ionisation_mode='pos', library_mode = 'standard', canopus_level = 'all', cosine=0.7, shared_peaks=6, ppm_error=10):
#table is a feature metadata table
#ionisation_mode = pos or neg
#library_mode = standard or analogue
#canopus_level = all or specific

    GNPS_SIRIUS_annotations = table
    
    #=============== ANALOGUE LIBRARY SEARCH PREPARATION ===========
    # For standard library search, we extra molecular formula from InChI and use it to evaluate ZODIAC results
    if library_mode.startswith("ana"):
        print('=== Looking at match between GNPS library in ANALOGUE mode and SIRIUS annotation ===')
        GNPS_SIRIUS_annotations = table
        #Keep only the entries with both GNPS INCHI available and a SIRIUS MF
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIBA_superclass.notnull() & 
                                                          GNPS_SIRIUS_annotations.SIR_MF_Zod_molecularFormula.notnull()]

        # First we filter likely incorrect spectral library match    
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIBA_IonMode.str.lower().str.startswith(ionisation_mode, na=False)] #check ionisation mode
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIBA_SpecCharge <= 1 ]
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIBA_SharedPeaks > shared_peaks]
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIBA_MQScore > cosine]

        usable_GNPS_SIRIUS = GNPS_SIRIUS_annotations.shape[0]
    
    #=============== ELSE WE LOOK AT STANDARD LIBRARY SEARCH ===========
    else:
        print('=== Looking at match between GNPS library in REGULAR mode and SIRIUS annotation ===')
        #Keep only the entries with both GNPS INCHI available and a SIRIUS MF
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_INCHI.notnull() & 
                                                      GNPS_SIRIUS_annotations.SIR_MF_Zod_molecularFormula.notnull()]
        total_GNPS_SIRIUS = GNPS_SIRIUS_annotations.shape[0]
        
        # First we filter likely incorrect spectral library match    
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_IonMode.str.lower().str.startswith(ionisation_mode, na=False)] #check ionisation mode
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_MZErrorPPM < ppm_error]
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_SpecCharge <= 1 ]
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_SharedPeaks > shared_peaks]
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_MQScore > cosine]

        # Entries that have salts or charged structure should be removed
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[~GNPS_SIRIUS_annotations.GNPS_LIB_INCHI.str.contains('q:+|p+1|p+2|p-1', na=False)]

        # Check that these are InChI
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_INCHI.str.upper().str.startswith(('INCHI','1S'), na=False)]

        # Filter the in source fragment spectral library match that are hard to match
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[~GNPS_SIRIUS_annotations.GNPS_LIB_Adduct.str.contains("-C|i")]

        # Get the MF from the INCHI string
        GNPS_SIRIUS_annotations['GNPS_LIB_INCHI_MF'] = GNPS_SIRIUS_annotations.GNPS_LIB_INCHI.str.split("/",expand=False).str[1]

        #Just to check there was no weird InChI
        GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.GNPS_LIB_INCHI_MF.str.upper().str.startswith('C', na=False)]
        usable_GNPS_SIRIUS = GNPS_SIRIUS_annotations.shape[0]
    
    # Number of entries with GNPS and SIRIUS
    print('Usable GNPS/SIRIUS annotations = '+str(usable_GNPS_SIRIUS))
    
    # Number of SIRIUS ZODIAC MOLECULAR FORMULA annotation passing the ZodiacScore
    GNPS_SIRIUS_annotations_score = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]
    print('Usable GNPS/SIRIUS annot. w. ZodiacScore > '+str(zodiac_score_thresh)+' = '+str(GNPS_SIRIUS_annotations.shape[0]))

    
    # ======= For regular library only. We check MF match
    if library_mode.startswith("ana") == False:
        # We are matching molecular formulas only by the number of carbons as the SIRIUS adducts can be quite variable and deals strangely with adduct/ionisationin the MF. Like loss of water etc. Ammonium adduct.
        GNPS_SIRIUS_annotations['GNPS_LIB_INCHI_MF_partial'] = GNPS_SIRIUS_annotations['GNPS_LIB_INCHI_MF'].str[:3]
        GNPS_SIRIUS_annotations['SIR_MF_Zod_molecularFormula_partial'] = GNPS_SIRIUS_annotations['SIR_MF_Zod_molecularFormula'].str[:3]

        # Make functions to check the match
        def mf_check(x):    
           return 'yes' if x['SIR_MF_Zod_molecularFormula_partial'] == x['GNPS_LIB_INCHI_MF_partial'] else 'no'

        # Apply the functions for all valid MF pairs above zodiac score
        GNPS_SIRIUS_annotations['MF_match'] = GNPS_SIRIUS_annotations.apply(mf_check, axis=1)
        MF_no_match = GNPS_SIRIUS_annotations[(GNPS_SIRIUS_annotations['MF_match'] == 'no')]
        MF_pairs = GNPS_SIRIUS_annotations.dropna(subset=['MF_match'])
        MF_match = GNPS_SIRIUS_annotations[(GNPS_SIRIUS_annotations['MF_match'] == 'yes')]
        
        # Apply the functions for all usable pairs above zodiac score
        GNPS_SIRIUS_annotations_score = GNPS_SIRIUS_annotations[GNPS_SIRIUS_annotations.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]        
        MF_no_match_score = GNPS_SIRIUS_annotations_score[(GNPS_SIRIUS_annotations_score['MF_match'] == 'no')]
        MF_pairs_score = GNPS_SIRIUS_annotations_score.dropna(subset=['MF_match'])
        MF_match_score = GNPS_SIRIUS_annotations_score[(GNPS_SIRIUS_annotations_score['MF_match'] == 'yes')]

                                                                                   
    # ==== CHEMICAL CLASSIFICATION CHECK
    # Checking if the GNPS lib class annotation in the SIRIUS classification
    def check_all_categories(table,column1substring,column2string,new_column_name):
        new_column = []
    
        for i, row in table.iterrows():
            if row[column1substring] is np.nan:
                new_column.append(np.nan)
                #print(str(row[column1substring])+' is np.nan')
            elif row[column2string] is np.nan:
                new_column.append(np.nan)
                #print(str(row[column2string])+ 'is np.nan')
            elif row[column1substring] not in row[column2string]:
                new_column.append('no')
                #print(str(row[column1substring])+ '=== IS NOT ===' + str(row[column2string]))
            elif row[column1substring] in row[column2string]:
                new_column.append('yes')
                #print(str(row[column1substring])+ '### IS ###' + str(row[column2string]))
            else:
                new_column.append(np.nan)
                #print(str(row[column1substring])+ '### OTHER ###' + str(row[column2string]))

        table[new_column_name] = new_column

        
    # ======= Check if classification is correct for ANALOGUE library search    
    if library_mode.startswith("ana"):
        gnps_lib = 'GNPS_LIBA_'
    else:
        gnps_lib = 'GNPS_LIB_'
        
    if canopus_level != 'all':
        print('Check with CANOPUS SPECIFIC classification levels')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'subclass'),'CAN_subclass',\
                                'Match_GNPSsubclass-SIRIUS')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'class'),'CAN_class',\
                                'Match_GNPSclass-SIRIUS')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'superclass'),'CAN_superclass',\
                                'Match_GNPSsuperclass-SIRIUS')     
    else:
        print('Check with CANOPUS ALL classification level')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'subclass'),'CAN_all classifications',\
                                'Match_GNPSsubclass-SIRIUS')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'class'),'CAN_all classifications',\
                                'Match_GNPSclass-SIRIUS')
        check_all_categories(GNPS_SIRIUS_annotations,(str(gnps_lib)+'superclass'),'CAN_all classifications',\
                                'Match_GNPSsuperclass-SIRIUS')

    print(' ')
    # Print MF results if regular library mode
    if library_mode.startswith('reg'):
        print('====== Match for molecular formulas =======')
        print('MF match = '+str(MF_match.shape[0]))
        print('MF match score = '+str(MF_match_score.shape[0]))
        print(' ')
        
    print('====== Match between GNPS lib superclass/class/subclass and SIRIUS CANOPUS level(s) =======')
    superclass_match_all_total = GNPS_SIRIUS_annotations.dropna(subset=['Match_GNPSsuperclass-SIRIUS'])
    superclass_match_all_total = superclass_match_all_total[superclass_match_all_total.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]
    superclass_match_all = superclass_match_all_total[(superclass_match_all_total['Match_GNPSsuperclass-SIRIUS'] == 'yes')]
    print('Classified pairs considered = '+str(superclass_match_all_total.shape[0]))
    print('Superclass annotation pairs = '+str(superclass_match_all.shape[0]))
    print('Superclass match all = '+str(superclass_match_all.shape[0])+', 'f"{superclass_match_all.shape[0]/superclass_match_all_total.shape[0]:.2f}"+'%')
    class_match_all_total = GNPS_SIRIUS_annotations.dropna(subset=['Match_GNPSclass-SIRIUS'])
    class_match_all_total = class_match_all_total[class_match_all_total.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]
    class_match_all = class_match_all_total[(class_match_all_total['Match_GNPSclass-SIRIUS'] == 'yes')]
    print('Class annotation pairs = '+str(class_match_all_total.shape[0]))
    print('Class match = '+str(class_match_all.shape[0])+', 'f"{class_match_all.shape[0]/class_match_all_total.shape[0]:.2f}"+'%')
    subclass_match_all_total = GNPS_SIRIUS_annotations.dropna(subset=['Match_GNPSsubclass-SIRIUS'])
    subclass_match_all_total = subclass_match_all_total[subclass_match_all_total.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]
    subclass_match_all = subclass_match_all_total[(subclass_match_all_total['Match_GNPSsubclass-SIRIUS'] == 'yes')]
    print('Subclass annotation pairs = '+str(subclass_match_all_total.shape[0]))
    print('Subclass match all = '+str(subclass_match_all.shape[0])+', 'f"{subclass_match_all.shape[0]/subclass_match_all_total.shape[0]:.2f}"+'%')
    
    check_matching_annotations.GNPS_SIRIUS_annotations = GNPS_SIRIUS_annotations

    check_matching_annotations.superclass_match_all_total = superclass_match_all_total
    check_matching_annotations.class_match_all_total = class_match_all_total
    check_matching_annotations.subclass_match_all_total = subclass_match_all_total

    check_matching_annotations.superclass_match_all = superclass_match_all
    check_matching_annotations.class_match_all = class_match_all
    check_matching_annotations.subclass_match_all = subclass_match_all


    #============
    if library_mode.startswith('reg'):
        data_matching = [['Usable MF pairs',MF_pairs.shape[0]],
                         ['Usable MF pairs w. ZodiacScore>'+str(zodiac_score_thresh), GNPS_SIRIUS_annotations_score.shape[0]],
                         ['Matching molecular formula', MF_match.shape[0]],
                         ['Matching molecular w. ZodiacScore>'+str(zodiac_score_thresh), MF_match_score.shape[0]]]
    
        check_matching_annotations.table_matching = pd.DataFrame(data_matching, columns = ['Matching level', 'Count'])
    
        #Make a column in relative values
        check_matching_annotations.table_matching['Relative'] = check_matching_annotations.table_matching['Count'].div(GNPS_SIRIUS_annotations_score.shape[0]).round(2)
        check_matching_annotations.MF_match = MF_match 
        check_matching_annotations.MF_no_match = MF_no_match
        check_matching_annotations.MF_pairs = MF_pairs
    
    else:
        data_matching = [['Usable pairs', usable_GNPS_SIRIUS],
                     ['Usable pairs w. ZodiacScore>'+str(zodiac_score_thresh), GNPS_SIRIUS_annotations_score.shape[0]]]
        
        check_matching_annotations.table_matching = pd.DataFrame(data_matching, columns = ['Matching level', 'Count'])
    
        #Make a column in relative values    
        check_matching_annotations.table_matching['Relative'] = check_matching_annotations.table_matching['Count'].div(usable_GNPS_SIRIUS).round(2)
            
    #=============
    #Make a table for CANOPUS
    data_class_matching = [['Available pairs', GNPS_SIRIUS_annotations.dropna(subset=['Match_GNPSsuperclass-SIRIUS']).shape[0]],
                     ['Classified pairs w. ZodiacScore>'+str(zodiac_score_thresh), superclass_match_all_total.shape[0]],
                     ['Matching superclass', superclass_match_all.shape[0]],
                     ['Matching class', class_match_all.shape[0]],['Matching subclass', subclass_match_all.shape[0]]]
    
    check_matching_annotations.table_class_matching = pd.DataFrame(data_class_matching, columns = ['Matching level', 'Count'])
    
    #Make a column in relative values    
    check_matching_annotations.table_class_matching['Relative'] = check_matching_annotations.table_class_matching['Count'].div(superclass_match_all_total.shape[0]).round(2)
 