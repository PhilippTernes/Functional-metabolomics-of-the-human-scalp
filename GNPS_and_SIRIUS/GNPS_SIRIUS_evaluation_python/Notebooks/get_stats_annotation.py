# Author Louis Felix Nothias Feb 2021
import numpy as np
import pandas as pd

# A function that gets stats for the GNPS annotations, Dereplicator, SIRIUS and also includes Passatutto
def get_stats_annotation(table, zodiac_score_thresh=0.7, ionisation_mode='pos', ppm_error=25):
    # ionisation_mode='pos' or 'neg'
    # Number of features
    df = table
    print('Features = '+str(df.shape[0]))

    print(' ')
    print('==== GNPS =====')
    # Number of features in molecular networks
    network = df[df.GNPS_componentindex.notnull()]
    network = network[network.GNPS_componentindex != -1]
    print('In networks = '+str(network.shape[0]))

    # Number of molecular networks
    number_of_networks = len(set(network['GNPS_componentindex']))
    print('Number of networks = '+str(number_of_networks))
    
    # Number of library annotation
    lib = df[df.GNPS_LIB_SpectrumID.notnull()]
    lib = lib[lib.GNPS_LIB_IonMode.str.lower().str.startswith(ionisation_mode, na=False)] #check ionisation mode
    lib = lib[lib.GNPS_LIB_MZErrorPPM < ppm_error]
    print('Valid library annotations = '+str(lib.shape[0]))

    # Number of library annotation in analogue mode
    liba = df[df.GNPS_LIBA_SpectrumID.notnull()]
    liba = liba[liba.GNPS_LIBA_IonMode.str.lower().str.startswith(ionisation_mode, na=False)] #check ionisation mode
    print('Library annotations in analogue mode= '+str(liba.shape[0]))

    # Number of passatutto annotation
    passa = df[df.PASSA_FDR_fdr.notnull()]
    passa = passa[passa.PASSA_FDR_IonMode.str.lower().str.startswith(ionisation_mode, na=False)] #check ionisation mode
    passa = passa[passa.PASSA_FDR_MZErrorPPM < ppm_error]
    print('PASSATUTTO FDR-controlled library annotations = '+str(passa.shape[0]))
    
    # Number of passatutto annotation FDR 0.8
    passa_score_10 = passa[passa.PASSA_FDR_fdr <= 0.1]
    passa_score_20 = passa[passa.PASSA_FDR_fdr <= 0.2]
    print('PASSATUTTO FDR-controlled library annotations at 20% FDR = '+str(passa_score_20.shape[0]))
    print('PASSATUTTO FDR-controlled library annotations at 10% FDR = '+str(passa_score_10.shape[0]))

    
    print(' ')
    print('==== SIRIUS =====')
    # Number of SIRIUS ZODIAC MOLECULAR FORMULA (MF) annotation
    sir_zod = df[df.SIR_MF_Zod_molecularFormula.notnull()]
    total_sir =  sir_zod.shape[0]
    print('Features with SIRIUS annotation = '+str(total_sir))

    # Number of SIRIUS ZODIAC MF annotation passing the ZodiacScore 
    sir_zod = sir_zod[sir_zod.SIR_MF_Zod_ZodiacScore > zodiac_score_thresh]
    print('SIRIUS ZODIAC MF with ZodiacScore > '+str(zodiac_score_thresh)+' = '+str(sir_zod.shape[0]))

    # Number of CSI FingerID annotations for MF passing the ZodiacScore 
    csi = sir_zod[sir_zod.CSI_smiles.notnull()]
    print('CSIFingerID annotations = '+str(csi.shape[0]))

    # Number of CANOPUS annotations for MF passing the ZodiacScore 
    can = sir_zod[sir_zod['CAN_all classifications'].notnull()]
    print('CANOPUS annotations = '+str(can.shape[0]))

    # Number of CANOPUS annotations up to the subclass level
    can_class = can[can['CAN_class'].notnull()]
    print('CANOPUS annotations at the subclass level= '+str(can_class.shape[0]))
    
    # Number of CANOPUS annotations up to the subclass level
    can_sub = can[can['CAN_subclass'].notnull()]
    print('CANOPUS annotations at the subclass level= '+str(can_sub.shape[0]))

    # Number of CANOPUS annotations up to the most specific level
    can_spe = can[can['CAN_level 5'].notnull()]
    print('CANOPUS annotations at the level 5 = '+str(can_spe.shape[0]))
    
    print(' ')
    print('==== General annotation statistics =====')

    #Number of features
    print('Number of features = '+ str(df.shape[0]))
    
    #Number of features with a least one annotation
    df_one_annotation = df.loc[(df.SIR_MF_Zod_molecularFormula.notnull())\
           |(df.GNPS_LIB_SpectrumID.notnull()) | (df.GNPS_LIBA_SpectrumID.notnull())]
    print('Annotated features = '+ str(df_one_annotation.shape[0]))
    
    #Number of features with a least one annotation or in network
    df_one_annotation_or_network = df.loc[((df.GNPS_componentindex != -1)) | (df.SIR_MF_Zod_molecularFormula.notnull())\
           |(df.GNPS_LIB_SpectrumID.notnull()) | (df.GNPS_LIBA_SpectrumID.notnull())\
           |(df.PASSA_FDR_fdr.notnull())]
    print('Annotated features or in network = '+ str(df_one_annotation_or_network.shape[0]))
    
    # Check annotation per networks    
    get_stats_annotation.network = network
    get_stats_annotation.network_dict_SIR = network.groupby('GNPS_componentindex')['SIR_MF_Zod_molecularFormula'].apply(list).to_dict()
    get_stats_annotation.network_dict_LIB = network.groupby('GNPS_componentindex')['GNPS_LIB_SpectrumID'].apply(list).to_dict()
    get_stats_annotation.network_dict_LIBA = network.groupby('GNPS_componentindex')['GNPS_LIBA_SpectrumID'].apply(list).to_dict()
    
    #['GNPS_LIB_SpectrumID','GNPS_LIBA_SpectrumID','DEREP_Score','DEREP+_Score']

    #Number of features that are single nodes
    df_singlenode = df.loc[(df.GNPS_componentindex == -1)]
    print('Single nodes = '+ str(df_singlenode.shape[0]))

    #Number of features that are single nodes and unannotated
    df_unannotated_and_singlenode = df.loc[(df.GNPS_componentindex == -1)
                                           & (df.GNPS_LIB_SpectrumID.isnull())
                                              & (df.GNPS_LIBA_SpectrumID.isnull())
                                              & (df.SIR_MF_Zod_molecularFormula.isnull())
                                              ]
    print('Single nodes and unnnannotated = '+ str(df_unannotated_and_singlenode.shape[0]))
    
    print(' ')
    
    #Make a big table with all informations
    total = df.shape[0]
    data = [['Features', total], ['GNPS - in networks', network.shape[0]], ['GNPS - lib. match', lib.shape[0]],
            ['GNPS - lib. match analogue', liba.shape[0]], 
            #['Unnannoted single nodes', df_unannotated_and_singlenode.shape[0]],
            ['PASSATUTTO FDR 20%', passa_score_20.shape[0]],['PASSATUTTO FDR 10%', passa_score_10.shape[0]],
            ['SIRIUS - Annotated features', total_sir], ['SIRIUS - MF with ZodScore >'+str(zodiac_score_thresh), sir_zod.shape[0]],
            ['SIRIUS - structure', csi.shape[0]],
            ['SIRIUS - chemical class', can.shape[0]], 
            ['Annotated features', df_one_annotation.shape[0]],
            ['Annotated features or in network', df_one_annotation_or_network.shape[0]]]
    get_stats_annotation.final_table = pd.DataFrame(data, columns = ['Annotation tool', 'Count'])
    get_stats_annotation.final_table_rel = pd.DataFrame(data, columns = ['Annotation tool', 'Count'])

    #Make GNPS table 
    data_gnps = [['Features', df.shape[0]], ['GNPS - in networks', network.shape[0]], ['GNPS - lib match', lib.shape[0]],
            ['GNPS - lib. match analogue ', liba.shape[0]], 
            ['FDR Passatutto 10%', passa_score_10.shape[0]],['FDR Passatutto 20%', passa_score_20.shape[0]]]
    get_stats_annotation.table_gnps = pd.DataFrame(data_gnps, columns = ['Annotation tool', 'Count'])
    get_stats_annotation.table_gnps_rel = pd.DataFrame(data_gnps, columns = ['Annotation tool', 'Count'])

    #Make a SIRIUS table
    data_sirius = [['Features', df.shape[0]], ['SIRIUS - Annotated features', total_sir],
            ['SIRIUS - MF at ZodScore >'+str(zodiac_score_thresh), sir_zod.shape[0]], ['SIRIUS - structure', csi.shape[0]],
            ['SIRIUS - superclass', can.shape[0]], ['SIRIUS - class', can_class.shape[0]],
            ['SIRIUS - subclass', can_sub.shape[0]]]
    get_stats_annotation.table_sirius = pd.DataFrame(data_sirius, columns = ['Annotation tool', 'Count'])
    get_stats_annotation.table_sirius_rel = pd.DataFrame(data_sirius, columns = ['Annotation tool', 'Count'])

    #Make a PASSATUTTO table
    data_passa = [['Features', df.shape[0]], ['FDR Passatutto 20%', passa_score_20.shape[0]],['FDR Passatutto 10%', passa_score_10.shape[0]]]
    get_stats_annotation.table_passa = pd.DataFrame(data_passa, columns = ['Annotation tool', 'Count'])
    get_stats_annotation.table_passa_rel = pd.DataFrame(data_passa, columns = ['Annotation tool', 'Count'])
    
    #Transform into proportions
    get_stats_annotation.final_table_rel['Count'] = get_stats_annotation.final_table_rel['Count'].div(total).round(2)
    get_stats_annotation.table_gnps_rel['Count'] =  get_stats_annotation.table_gnps_rel['Count'].div(total).round(2)
    get_stats_annotation.table_sirius_rel['Count'] = get_stats_annotation.table_sirius_rel['Count'].div(total).round(2)
    get_stats_annotation.table_passa_rel['Count'] = get_stats_annotation.table_passa_rel['Count'].div(total).round(2)