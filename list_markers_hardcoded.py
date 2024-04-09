




nsforest0_thalassa = ['Atp6v0d2', 'Abcg1',# AM
             'Rtkn2',  'Igfbp2', #AT1
             'Sftpc','Cxcl15', #AT2,
            'Cd79a', #B_cells
             'Ms4a2', 'Fcer1a', #Basophils
             'Ccdc153', #Ciliated
             'Scgb3a2', 'Scgb1a1',#Club
             'Cst3',#DC
             'Cdh5', 'Clec14a',  #EC
             'Inmt', 'Pcolce2', # Fibroblasts
             'C1qc', 'C1qa', 'C1qb', # 'C3ar1', #IM
             'Upk3b',# Mesotheliocytes
             'Ifitm6','Plac8',# Monocytes
            'Ms4a4b', 'Ccl5', 'Hcst', # NK_T_cells
             'Gzma', 'Ncr1',# NK_cells
             'S100a9',# Neutrophils
             'Mmrn1',#Platelets
           'Acta2','Myh11', # SMC
             'Cd3g', 'Cd3d' #T_cells
             ]



list_cube0  = ['A', "B"]
list_cubeAB  = ['A', "B"]

list_cubeABC  = ['A', "B", "C"]
list_cubeABCD  = ['A', "B", "C", "D"]



#### color code


palette = {
    "-1": "black",
    "0": "blue",
    "1": "orange",
    "2": "green",
    "3": "red",
    "4": "purple",
    "5": "brown",
    "6": "pink",
    "7": "olive",
    "8": "cyan",
    "9": "darkseagreen",
    "10": "coral",
    "11": 'lime',
    "12": "gold",
    "13": "lime",
    "14": "darkseagreen",
    "15": "chocolate",
    "16": "coral", "17": "plum", "18": "magenta",
    "19": "rosybrown", "20": "khaki", "21": "teal"
}



oligo_pool2_6g = ['Rtkn2',
                 'Pecam1',
                  'Ptprb',
                  'Pdgfra',
                  'Apln',
                  'Hhip']





iss_topo141  = ['ARIH1', 'ATP11A', 'BCL2', 'BMP5', 'CCBE1', 'CCL21', 'CD36',
       'CD74', 'CD93', 'CDON', 'CLDN5', 'COL11A1', 'COL13A1', 'COL9A1',
       'COL9A3', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1', 'DLL3', 'DLL4',
       'DMD', 'DNAH12', 'DTX1', 'EBF1', 'EDN1', 'EGFL6', 'EPCAM', 'ETS1',
       'ETS2', 'ETV1', 'ETV3', 'ETV5', 'ETV6', 'FAM162B', 'FGF10',
       'FGF18', 'FGF2', 'FGF20', 'FGF7', 'FGF9', 'FGFR1', 'FGFR2',
       'FGFR3', 'FGFR4', 'FLT1', 'FLT4', 'FZD1', 'FZD2', 'FZD3', 'FZD6',
       'FZD7', 'G0S2', 'GLI2', 'GLI3', 'GPC5', 'GRP', 'HES1', 'HEY1',
       'HEYL', 'HGF', 'HHIP', 'IGF1', 'IGFBP5', 'IGFBP7', 'ITGA8', 'JAG1',
       'JAG2', 'KCNQ5', 'KDR', 'KLRB1', 'LEF1', 'LRP1B', 'LRP2', 'LRP5',
       'LTB', 'LXN', 'LYZ', 'MEOX2', 'MET', 'MFNG', 'MGP', 'MYH11',
       'NKD1', 'NKG7', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOTUM',
       'PDGFA', 'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PDZRN4',
       'PECAM1', 'PHOX2B', 'PRKCB', 'PTCH1', 'RELN', 'RGS7', 'RNASE1',
       'RSPH1', 'RSPO2', 'RSPO3', 'S100A9', 'SAMD4A', 'SEC11C', 'SFRP1',
       'SFRP2', 'SFTPC', 'SHH', 'SMO', 'SORCS1', 'SOX17', 'SPOPL',
       'SPRY1', 'SPRY2', 'SRGN', 'STAT3', 'STMN2', 'TAGLN', 'TCF7L1',
       'TCIM', 'TECRL', 'TNNT1', 'TPPP3', 'UCMA', 'UPK3B', 'VEGFC',
       'VEGFD', 'WIF1', 'WNT11', 'WNT2', 'WNT2B', 'WNT5A', 'WNT7B', 'WT1']


iss_topo142 = ['AP000561', 'ARIH1', 'ATP11A', 'BCL2', 'BMP5', 'CCBE1', 'CCL21', ## for pwc6
       'CD36', 'CD74', 'CD93', 'CDON', 'CLDN5', 'COL11A1', 'COL13A1',
       'COL9A1', 'COL9A3', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1',
       'DLL3', 'DLL4', 'DMD', 'DNAH12', 'DTX1', 'EBF1', 'EDN1', 'EGFL6',
       'EPCAM', 'ETS1', 'ETS2', 'ETV1', 'ETV3', 'ETV5', 'ETV6', 'FAM162B',
       'FGF10', 'FGF18', 'FGF2', 'FGF20', 'FGF7', 'FGF9', 'FGFR1',
       'FGFR2', 'FGFR3', 'FGFR4', 'FLT1', 'FLT4', 'FZD1', 'FZD2', 'FZD3',
       'FZD6', 'FZD7', 'G0S2', 'GLI2', 'GLI3', 'GPC5', 'GRP', 'HES1',
       'HEY1', 'HEYL', 'HGF', 'HHIP', 'IGF1', 'IGFBP5', 'IGFBP7', 'ITGA8',
       'JAG1', 'JAG2', 'KCNQ5', 'KDR', 'KLRB1', 'LEF1', 'LRP1B', 'LRP2',
       'LRP5', 'LTB', 'LXN', 'LYZ', 'MEOX2', 'MET', 'MFNG', 'MGP',
       'MYH11', 'NKD1', 'NKG7', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4',
       'NOTUM', 'PDGFA', 'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB',
       'PDZRN4', 'PECAM1', 'PHOX2B', 'PRKCB', 'PTCH1', 'RELN', 'RGS7',
       'RNASE1', 'RSPH1', 'RSPO2', 'RSPO3', 'S100A9', 'SAMD4A', 'SEC11C',
       'SFRP1', 'SFRP2', 'SFTPC', 'SHH', 'SMO', 'SORCS1', 'SOX17',
       'SPOPL', 'SPRY1', 'SPRY2', 'SRGN', 'STAT3', 'STMN2', 'TAGLN',
       'TCF7L1', 'TCIM', 'TECRL', 'TNNT1', 'TPPP3', 'UCMA', 'UPK3B',
       'VEGFC', 'VEGFD', 'WIF1', 'WNT11', 'WNT2', 'WNT2B', 'WNT5A',
       'WNT7B', 'WT1']


iss_topo130_pcw13 = ['AP000561', 'ARIH1', 'ATP11A', 'BCL2', 'BMP5', 'CCBE1', 'CCL21',
       'CD36', 'CD74', 'CD93', 'CLDN5', 'COL11A1', 'COL13A1', 'COL9A1',
       'COL9A3', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1', 'DLL3', 'DLL4',
       'DMD', 'DNAH12', 'DTX1', 'EBF1', 'EDN1', 'EGFL6', 'EPCAM', 'ETS1',
       'ETS2', 'ETV1', 'ETV3', 'ETV5', 'ETV6', 'FAM162B', 'FGF10',
       'FGF18', 'FGF2', 'FGF20', 'FGF7', 'FGFR1', 'FGFR2', 'FGFR3',
       'FGFR4', 'FLT4', 'FZD1', 'FZD2', 'FZD3', 'FZD6', 'FZD7', 'G0S2',
       'GLI2', 'GLI3', 'GPC5', 'GRP', 'HES1', 'HEY1', 'HEYL', 'HGF',
       'HHIP', 'IGF1', 'IGFBP5', 'ITGA8', 'JAG1', 'JAG2', 'KCNQ5', 'KDR',
       'KLRB1', 'LEF1', 'LRP1B', 'LRP2', 'LRP5', 'LTB', 'LXN', 'MEOX2',
       'MET', 'MFNG', 'MGP', 'MYH11', 'MYOCD', 'NKD1', 'NKG7', 'NOTCH1',
       'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOTUM', 'NRXN1', 'PDGFA', 'PDGFB',
       'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PDZRN4', 'PHOX2B', 'PRKCB',
       'PTCH1', 'RELN', 'RGS7', 'RNASE1', 'RSPH1', 'RSPO2', 'S100A9',
       'SAMD4A', 'SHH', 'SMO', 'SORCS1', 'SOX17', 'SPOPL', 'SPRY1',
       'SPRY2', 'SRGN', 'STAT3', 'STMN2', 'TCF7L1', 'TCIM', 'TECRL',
       'TNNT1', 'TPPP3', 'UCMA', 'VEGFB', 'VEGFC', 'WIF1', 'WNT11',
       'WNT2', 'WNT5A', 'WNT7B', 'WT1']


iss_topo129_pcw13 = ['AP000561', 'ARIH1', 'ATP11A', 'BCL2', 'BMP5', 'CCBE1', 'CCL21',
       'CD36', 'CD74', 'CD93', 'CLDN5', 'COL11A1', 'COL13A1', 'COL9A1',
       'COL9A3', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1', 'DLL3', 'DLL4',
       'DMD', 'DNAH12', 'DTX1', 'EDN1', 'EGFL6', 'EPCAM', 'ETS1', 'ETS2',
       'ETV1', 'ETV3', 'ETV5', 'ETV6', 'FAM162B', 'FGF10', 'FGF18',
       'FGF2', 'FGF20', 'FGF7', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4',
       'FLT4', 'FZD1', 'FZD2', 'FZD3', 'FZD6', 'FZD7', 'G0S2', 'GLI2',
       'GLI3', 'GPC5', 'GRP', 'HES1', 'HEY1', 'HEYL', 'HGF', 'HHIP',
       'IGF1', 'IGFBP5', 'ITGA8', 'JAG1', 'JAG2', 'KCNQ5', 'KDR', 'KLRB1',
       'LEF1', 'LRP1B', 'LRP2', 'LRP5', 'LTB', 'LXN', 'MEOX2', 'MET',
       'MFNG', 'MGP', 'MYH11', 'MYOCD', 'NKD1', 'NKG7', 'NOTCH1',
       'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOTUM', 'NRXN1', 'PDGFA', 'PDGFB',
       'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PDZRN4', 'PHOX2B', 'PRKCB',
       'PTCH1', 'RELN', 'RGS7', 'RNASE1', 'RSPH1', 'RSPO2', 'S100A9',
       'SAMD4A', 'SHH', 'SMO', 'SORCS1', 'SOX17', 'SPOPL', 'SPRY1',
       'SPRY2', 'SRGN', 'STAT3', 'STMN2', 'TCF7L1', 'TCIM', 'TECRL',
       'TNNT1', 'TPPP3', 'UCMA', 'VEGFB', 'VEGFC', 'WIF1', 'WNT11',
       'WNT2', 'WNT5A', 'WNT7B', 'WT1']


iss_topo89_pcw13 = ['BCL2', 'CCL21', 'CDON', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1', 'DLL3', 'DLL4', 'DNAH12',
                    'DTX1', 'ETS1', 'ETS2', 'ETV1', 'ETV3', 'ETV5', 'ETV6',
                    'FGF10', 'FGF18', 'FGF2', 'FGF20', 'FGF7', 'FGF9', 'FGFR1', 'FGFR2', 'FGFR3', 'FGFR4', 'FLT1',
                    'FLT4', 'FZD1', 'FZD2', 'FZD3', 'FZD6', 'FZD7', 'GLI2',
                    'GLI3', 'GRP', 'HES1', 'HEY1', 'HEYL', 'HGF', 'IGFBP7', 'JAG1', 'JAG2', 'KDR', 'LEF1', 'LRP2',
                    'LRP5', 'LYZ', 'MET', 'MFNG', 'NKD1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOTUM', 'PDGFA',
                    'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PHOX2B', 'PTCH1', 'RSPH1', 'RSPO2', 'RSPO3', 'SEC11C', 'SFRP1',
                    'SFRP2', 'SHH', 'SMO', 'SOX17', 'SPOPL', 'SPRY1', 'SPRY2', 'TCF7L1', 'TPPP3', 'VEGFB', 'VEGFC',
                    'VEGFD', 'WIF1', 'WNT11', 'WNT2', 'WNT2B', 'WNT5A', 'WNT7B']



iss_topo138_pcw13 = ['AP000561', 'ARIH1', 'ATP11A', 'BCL2', 'BMP5', 'CCBE1', 'CCL21',
       'CD36', 'CD74', 'CD93', 'CLDN5', 'COL11A1', 'COL13A1', 'COL9A1',
       'COL9A3', 'CPM', 'CTGF', 'CTNND2', 'CXXC4', 'DLL1', 'DLL3', 'DLL4',
       'DMD', 'DNAH12', 'DTX1', 'EBF1', 'EDN1', 'EGFL6', 'EPCAM', 'ETS1',
       'ETS2', 'ETV1', 'ETV3', 'ETV5', 'ETV6', 'FAM162B', 'FGF10',
       'FGF18', 'FGF2', 'FGF20', 'FGF7', 'FGFR1', 'FGFR2', 'FGFR3',
       'FGFR4', 'FLT1', 'FLT4', 'FZD1', 'FZD2', 'FZD3', 'FZD6', 'FZD7',
       'G0S2', 'GLI2', 'GLI3', 'GPC5', 'GRP', 'HES1', 'HEY1', 'HEYL',
       'HGF', 'HHIP', 'IGF1', 'IGFBP5', 'ITGA8', 'JAG1', 'JAG2', 'KCNQ5',
       'KDR', 'KLRB1', 'LEF1', 'LRP1B', 'LRP2', 'LRP5', 'LTB', 'LXN',
       'MEOX2', 'MET', 'MFNG', 'MGP', 'MYH11', 'MYOCD', 'NKD1', 'NKG7',
       'NKX2-1', 'NOTCH1', 'NOTCH2', 'NOTCH3', 'NOTCH4', 'NOTUM', 'NRXN1',
       'PDGFA', 'PDGFB', 'PDGFC', 'PDGFD', 'PDGFRA', 'PDGFRB', 'PDZRN4',
       'PECAM1', 'PHOX2B', 'PRKCB', 'PTCH1', 'RELN', 'RGS7', 'RNASE1',
       'RSPH1', 'RSPO2', 'S100A9', 'SAMD4A', 'SFRP1', 'SFTPC', 'SHH',
       'SMO', 'SORCS1', 'SOX17', 'SPOPL', 'SPRY1', 'SPRY2', 'SRGN',
       'STAT3', 'STMN2', 'TAGLN', 'TCF7L1', 'TCIM', 'TECRL', 'TNNT1',
       'TPPP3', 'UCMA', 'UPK3B', 'VEGFB', 'VEGFC', 'VEGFD', 'WIF1',
       'WNT11', 'WNT2', 'WNT5A', 'WNT7B', 'WT1']


iss_topo81_pcw13 = ['PHOX2B', 'DLL4', 'FLT4', 'NOTCH1', 'DLL1', 'PDGFB', 'LRP2', 'NOTCH3', 'RSPH1', 'TCF7L1', 'SPRY1', 'CTNND2', 'FZD2', 'NOTUM', 'TPPP3',
 'ETV1', 'PDGFD', 'WNT5A', 'WIF1', 'VEGFB', 'KDR', 'GLI2', 'PDGFRB', 'FGF20', 'NOTCH4', 'VEGFC', 'LEF1', 'MFNG', 'FGFR3', 'PDGFRA', 'HGF',
 'SPOPL', 'FGFR4', 'DTX1', 'NKD1', 'ETV6', 'BCL2', 'SPRY2', 'CXXC4', 'PDGFA', 'FZD6', 'HEY1', 'ETS2', 'RSPO2', 'FZD7', 'GLI3', 'FZD3', 'FGFR2',
 'ETV5', 'CTGF', 'CPM', 'FGF2', 'MET', 'ETS1', 'PTCH1', 'WNT2', 'FLT1', 'WNT7B', 'JAG2', 'SHH', 'FGF18', 'NOTCH2', 'SMO', 'PDGFC', 'DLL3', 'FGF10',
 'FZD1', 'WNT11', 'DNAH12', 'HEYL', 'SFRP1', 'VEGFD', 'CCL21', 'LRP5', 'FGF7', 'FGFR1', 'GRP', 'ETV3', 'JAG1', 'SOX17', 'HES1']



list_cmap1 = ['k', 'k', '#fc8d59', '#ffffbf', '#91bfdb', '#19F7C1', '#CB9540', '#7938E5', '#A54AD2', '#EB9CF8', '#7BAA95', '#7C28BB', '#876DBE', '#1D9A63', '#5BC24E', '#4C3A22', '#85001D', '#5D9B6A', '#325675', '#FECD91', '#C90517', '#B5E55A', '#56D2EE', '#E14BB0', '#84EE0C', '#6FA54F', '#D276DB', '#2FBA17', '#CA5532', '#8F5B69', '#716863', '#323DA1', '#5E2527', '#41E821', '#779084', '#8E1EE6', '#F0E332', '#E96542', '#242985', '#5031D8', '#C93562', '#4F631D', '#0B55A4', '#E9F724', '#27D57E', '#353F72', '#E1A474', '#CEA6BC', '#4BD7B4', '#3FB674', '#E3E9A4', '#B5C7BC', '#41036A', '#AE3C28', '#7BC151', '#01C057', '#BBBA8C', '#4964D6', '#A0AB4C', '#4C9BD7', '#495BD3', '#36681A', '#AFDD94', '#E76CBA', '#E87047', '#5C773F', '#9C0F0B', '#F83561', '#F366CC', '#FFCAEF', '#F390C6', '#CBB8C2', '#8F0143', '#194378', '#4129A9', '#6871CD', '#124C28', '#844CA2', '#B59BE8', '#E78616', '#8F56CE', '#9AB2D9', '#FD7B76', '#3D1E87', '#90D56C', '#AC49D7', '#E243B4', '#436688', '#36D55B', '#AAC4A2', '#31F5B6', '#63EC38', '#55AEF8', '#4DC254', '#9123C3', '#777C50', '#50B5D0', '#152926', '#58A86E', '#C8E455', '#3568CA', '#8B1AB6', '#DBB14B', '#695A6C', '#E647B9', '#7FB45D', '#C1644E', '#69EBE0', '#7917D2', '#A8AC77', '#4F969C', '#B145F7', '#872DFD', '#B351C9', '#A3B8FC', '#24ECE8', '#D5E725', '#EDB361', '#0D7E0C', '#BEA222', '#89E86F', '#8C60B9', '#FD9C9A', '#12D74D', '#3F0971', '#9A3116', '#CE8531', '#8AEA91', '#FF333B', '#846B78', '#087ED6', '#BF27B6', '#737DE6', '#B7133E', '#382922', '#D572E5', '#476954', '#AA9403', '#1BAAA5', '#C00A5B', '#042752', '#FB517B', '#5D423B', '#D75E48', '#4C97E0', '#A4E5D5', '#508B57', '#8A96C8', '#F7EF55', '#5634B3', '#963876', '#2D0E2A', '#8426FB', '#40452B', '#222E7A', '#04B5F5', '#96035C', '#742FD5', '#331328', '#ED4EFF', '#663150', '#812219', '#557201', '#BAF73B', '#5AE30E', '#F97BB4', '#3BB27B', '#7F17E2', '#E90EBE', '#ED5E55', '#31176F', '#A7BA9A', '#1EF17E', '#051777', '#423725', '#BBE6F5', '#B7C42F', '#9AEF10', '#629B3C', '#70AFBA', '#A70ECE', '#646BC6', '#D30E5A', '#796AE1', '#DE4B90', '#DEDC19', '#1FB3A7', '#B1A733', '#F3F8C6', '#0C72F4', '#BBF8B2', '#B62887', '#6AFF4B', '#A7620B', '#8DDDB4', '#1E253D', '#EB5EAA', '#4B7D4D', '#E80006', '#4520BC', '#4B89F1', '#6569A4', '#08315A', '#7CBB65', '#C2D85E', '#AB643A', '#65915A', '#3DBB66', '#3C9EE9', '#A628E4', '#5F945C', '#206538', '#EA0543', '#CC01F4', '#5474A5', '#6492B8', '#586A99', '#19ABE7', '#B9213E', '#A08BE8', '#DA1775', '#01E57D', '#7914B8', '#1B7793', '#8FEE5F', '#5B25B1', '#4EB28F', '#55B11F', '#806F6A', '#FAD5CF', '#E1396E', '#D82DBC', '#E1895C', '#2DA174', '#57FE22', '#6278E5', '#47F649', '#31C5BF', '#C6E268', '#2C0092', '#AEF112', '#E2EAB4', '#7E4E45', '#F60E89', '#4B0400', '#EE2434', '#9178F6', '#CA5F05', '#9DCA6D', '#3C7BEF', '#82F57D', '#B67102', '#1A71F5', '#FC7B50', '#1563DC', '#61DFC6', '#A3DAAA', '#7B683E', '#C9628E', '#74F3A1', '#DA417E', '#303BF5', '#7B00F8', '#8EBB75', '#BAB67F', '#855D69', '#1443F1', '#DD070F', '#E1DBD3', '#171AC1', '#8496AF', '#6062FC', '#73076D', '#B7AB1C', '#FFE000', '#36DCA7', '#503913', '#96E785', '#DE4A1F', '#AB4B1B', '#3D31CE', '#8D69D6', '#D45A91', '#2DC00F', '#7B6BA6', '#073B10', '#F2D22B', '#DEAE10', '#527E48', '#533F43', '#D627BA', '#7BCCA1', '#7E3893', '#D0B146', '#E4E00D', '#8679D8', '#A5FC04', '#A1040B', '#7509ED', '#64D431', '#7A9B94', '#6096E3', '#88D269', '#DBBA14', '#75F3AC', '#A5080C', '#81B418', '#9B3F1F', '#32AA53', '#164FE2', '#EFB8E5', '#EF43F4', '#07EC0B', '#4E9E1F', '#F1A874', '#B92A8E', '#C285CB', '#052FFC', '#8772D8', '#104812', '#418D64', '#F0363B', '#802695', '#CBA0F7', '#68D663', '#5FF10D', '#46CC20', '#B99781', '#7A7C7D', '#D0FA14', '#3CFEA3', '#30B2F0', '#C5F7DD', '#1F37F5', '#EE4B29', '#193057', '#654767', '#E8DA60', '#8537F6', '#60065E', '#235EDD', '#18E40F', '#576AC2', '#D1800E', '#428AAF', '#86FEEE', '#D46883', '#A753CA', '#BD9054', '#288D75', '#56DBA6', '#BA38FA', '#DFF646', '#B3E708', '#478998', '#8B126B', '#426ED5', '#34DA7D', '#B7CFA5', '#E7BFA0', '#79A6AF', '#0914FC', '#3A2DA0', '#5D1EFE', '#E24432', '#82E835', '#E445A6', '#B47F1E', '#9E0A49', '#32637B', '#8EF1AF', '#34F6A6', '#D5350F', '#04510C', '#D6EF84', '#A8CD54', '#DBA90C', '#7426AE', '#C79BC2', '#FE9A28', '#A17232', '#5F0954', '#0919DF', '#C50074', '#14E21E', '#902BC1', '#68D020', '#F40394', '#2D83D2', '#EF604B', '#63AA3A', '#FC279F', '#1E3DC5', '#565170', '#25DC82', '#8F5932', '#EE8262', '#A537BE', '#4780DF', '#A045FF', '#5093F7', '#D11CFA', '#C4D5C1', '#494C9D', '#5D36F1', '#25CA4B', '#E46202', '#EBE537', '#00738A', '#67ADE3', '#DDF996', '#CA84FB', '#F1743D', '#6FFDAA', '#031334', '#D6D133', '#737918', '#9F4FB2', '#9BAC56', '#2C07DE', '#B01C97', '#037990', '#C55307', '#4BF891', '#C03971', '#4ACC0B', '#0F8286', '#CF9BDB', '#4D70DE', '#11B30D', '#06691F', '#A86A57', '#0E9D3D', '#7F96F3', '#89314E', '#EACFB9', '#40463A', '#0F8C24', '#36E7B6', '#4EE67B', '#C6150F', '#376F61', '#DF2D82', '#3A5869', '#414DB5', '#FD135D', '#A212F5', '#FD9D78', '#A40C12', '#0AFA7D', '#5140D0', '#29C6D5', '#9852D8', '#655B3B', '#F50FF3', '#AAC31A', '#D498FB', '#F3F362', '#8E51ED', '#E8A250', '#114692', '#B65A62', '#3E1F9A', '#7F765D', '#BB24D2', '#40F9A4', '#9B3DB7', '#71709D', '#325DAE', '#AFA82C', '#4CEECB', '#2BD334', '#BFA4CF', '#8FC188', '#6C11A3', '#6FCDE1', '#707C40', '#7AF7C9', '#31C181', '#C465B4', '#A4C5EC', '#6AAC1D', '#BC8781', '#6B9761', '#85A398', '#B3E720', '#B5D155', '#5DC6A3', '#9AE5B9', '#C07F0B', '#507233', '#5DC5AD', '#2C84D7', '#AAF015', '#2226CB', '#520257', '#03FEB4', '#68BB77', '#4082DB', '#C9C233', '#F86D24', '#894035', '#98B41C', '#9C1D7C', '#198138', '#FEB362', '#3CF994', '#830F1F', '#179F8F', '#E65DDC', '#E43ABE', '#34FF26', '#FD8969', '#C7247E', '#120B15', '#209520', '#67512A', '#A64586', '#0607C4', '#8DA1E3', '#DA3BBC', '#1FD10A', '#85C8D1', '#8CBBF4', '#DF0272', '#1E2436', '#31C016', '#42F274', '#5F6F32', '#E4C995', '#602269', '#AA5AA3', '#575868', '#01E077', '#3F539E', '#75D318', '#44407A', '#5384E7', '#D8A912', '#9D1CCF', '#980E03', '#343692', '#2BFC56', '#DECD85', '#15EF1C', '#82844F', '#606E86', '#6E3BFA', '#1A25BA', '#BEC00F', '#C48AA3', '#A81D03', '#6FFEE1', '#961229', '#D27F12', '#36A290', '#67C10B', '#0E742D', '#DFECF0', '#3FA2C0', '#17976E', '#860729', '#DACC3D', '#A3E015', '#3D8A72', '#A62A5D', '#CA9FFF', '#323FEE', '#BABC29', '#CEED81', '#FBD0B1', '#0915CF', '#717066', '#FCA174', '#A513B8', '#053FEE', '#32FBF0', '#848F96', '#A6902A', '#CAA5FA', '#2D7C20', '#8C7832', '#BAF7B5', '#64A1E9', '#F4C889', '#7805CD', '#82DA53', '#2E4C62', '#5D23E6', '#BAF3C0', '#BF596D', '#583B88', '#880EE6', '#B66FBC', '#2BF902', '#05CE3C', '#AA3A6E', '#F8894A', '#85475C', '#521E70', '#A2A994', '#0FAA6E', '#AFFD5C', '#DB97B6', '#6D8A59', '#426443', '#DEA07F', '#824746', '#D89348', '#DB8455', '#C6F4F2', '#47599E', '#104A3B', '#6FBC9D', '#3B3914', '#158A65', '#2356EF', '#05ABF3', '#34D6D5', '#4A2875', '#45C3BC', '#2041C0', '#DF7C84', '#A9EDC6', '#603206', '#38BD75', '#313957', '#5CA295', '#FCACC8', '#607153', '#DA32A6', '#2409F8', '#79A379', '#FE4F88', '#0A65D0', '#8B6A25', '#96A12F', '#8D4A35', '#50A472', '#41E009', '#B41937', '#453A45', '#1101E2', '#815212', '#5DFB26', '#2DE430', '#4686A1', '#A094AF', '#D45849', '#3AF541', '#050847', '#5A81F9', '#868C75', '#8647C5', '#DF9856', '#391071', '#7424AC', '#DF05B7', '#1DD5B1', '#6FCAE4', '#4A40E7', '#8876E8', '#08E5C1', '#ECCF46', '#AB0EA9', '#D048A5', '#DCE3A8', '#2B3770', '#DC79CE', '#2D0D85', '#5BF989', '#EA8E3E', '#83A630', '#C132F4', '#4D47DA', '#09BDFB', '#742473', '#7EB8AF', '#E2571D', '#03F443', '#10F543', '#3F583E', '#0B024E', '#4D9C67', '#0A63F9', '#2CE239', '#51169E', '#DEF212', '#076EFA', '#620B82', '#63F47F', '#DD16CE', '#D62053', '#D47248', '#A40750', '#25A443', '#4F1D51', '#E23DE4', '#3A1D79', '#7AEBE2', '#430AEC', '#7299CD', '#CF530E', '#56F595', '#86DBDD', '#3783AA', '#935CE6', '#26C1B6', '#6EDD9B', '#ECAB9A', '#3587AA', '#3B14AD', '#FE3142', '#F7A17B', '#899A7A', '#823E6B', '#6DCCF0', '#01ADB9', '#A524AA', '#3E0DDD', '#53C4E3', '#50F799', '#7D653B', '#20941B', '#A0B1C6', '#30B550', '#D471BF', '#FB4C5C', '#354B1D', '#68B888', '#08EB5A', '#31BB30', '#B9FEFB', '#F38E53', '#27780C', '#4BC0C6', '#90D54A', '#1229A0', '#F155C4', '#F11BCA', '#37CBD1', '#71AD03', '#19D626', '#C6F08D', '#26FF87', '#EEFB90', '#14A7DE', '#4D4431', '#1DDF64', '#EDA185', '#A75D38', '#1E6ECE', '#E9AF36', '#53166D', '#902E48', '#CAFB59', '#E47F86', '#B8D9BE', '#AAE814', '#4DBB34', '#3E9EAC', '#DC894D', '#F6868A', '#F781E1', '#4FC437', '#DC9186', '#CF6F1F', '#8A8F7C', '#CC640E', '#E7843E', '#E930A0', '#95D41A', '#280FDE', '#F931A2', '#326D2B', '#6144CC', '#E415FB', '#C62CB9', '#34DF0B', '#D472C7', '#380C25', '#850B06', '#5EDEA8', '#F708CB', '#8CE267', '#2BC0D1', '#86999A', '#A32679', '#C189FB', '#A3CA78', '#B0B1DA', '#0C4107', '#E02906', '#B2F650', '#A661CA', '#AC6A31', '#61F075', '#85EA77', '#6FA800', '#8404EE', '#7ADA83', '#1CC5E6', '#AE312B', '#30A450', '#BF33D8', '#E98644', '#4FECE8', '#73E7C4', '#C32255', '#932E9E', '#AE1FDD', '#E37080', '#D99918', '#C7DEA1', '#481D8D', '#6FA0D3', '#CB2955', '#ACC372', '#AA8BFA', '#30FF3B', '#43481F', '#C67D3D', '#F6F55A', '#0FFD84', '#794453', '#A61349', '#255D15', '#C91C63', '#419C0A', '#28E58B', '#43261A', '#07D4D4', '#9F62E0', '#2833AE', '#12CF60', '#ED2BD8', '#D50218', '#C75A0F', '#71DF3B', '#0E97C5', '#2BAD31', '#102E99', '#72257E', '#8230FB', '#6C8931', '#F05709', '#9DC9CA', '#749869', '#8AA780', '#5AD336', '#3121ED', '#8B84BF', '#4A07F5', '#2B8AA1', '#849684', '#B2D62D', '#C9FCFE', '#C9D367', '#CD9618', '#CD7143', '#3FC699', '#512E6E', '#A93E3A', '#98B532', '#8E3118', '#25B220', '#F1EB2D', '#FA0AE2', '#0A9D97', '#E67AE8', '#61151C', '#9B4CAF', '#4198C8', '#A20546', '#A0FB6A', '#F643F3', '#D31BE9', '#47CEE9', '#BA3AE0', '#192C17', '#40C533', '#DA7348', '#AFFDCF', '#4D403D', '#33B672', '#8C144B', '#023FE1', '#73F7A9', '#17B9C7', '#C75A61', '#BFCF61', '#15D6ED', '#1D9923', '#FF10E3', '#A7CA00', '#8C4579', '#EBAFBD', '#7CC7FD', '#1E03BD', '#038089', '#8C36E8', '#E3DE0D', '#5F7A7E', '#7178B5', '#0289B5', '#089FD7', '#160090', '#A35BFF', '#434E19', '#9AA4CA', '#4632FB', '#456CB2', '#204462', '#22B1A4', '#F6DF11', '#FF0101', '#EF660A', '#1AA0C5', '#8F7DE2', '#9A3725', '#6D823A', '#078DE9', '#91559D', '#F1FEE6', '#DB6F59', '#BDE37D', '#2F4A0F', '#52D509', '#13C8DD', '#C52D1C', '#54EEE6', '#9E1D21', '#901726', '#3C741C', '#EAF09B', '#B67204', '#A1FCAD', '#D0C072', '#327A5F', '#01B85D', '#0EB990', '#263995', '#08E677', '#689E92', '#804E72', '#F2BF78', '#2F2B80', '#BA3017', '#3C70A4', '#5434EF', '#A1F7AC', '#3A5754', '#CD235F', '#6C5139', '#201EB9', '#C5B34E', '#3E06A0', '#435116', '#17DD7A', '#4D1E7E', '#C57144', '#7A838B', '#C15530', '#C48FCB', '#BBD41B', '#08F6B4', '#FA1D40', '#6194FB', '#D5110A', '#5EAEA9', '#B4CB8A', '#50DDDC', '#C1A0E5', '#2FA2D0', '#0CCF9B', '#55B0F3', '#AB5436', '#C9ED19', '#12F802', '#50D958', '#31E8F1', '#CF3625', '#F70A64', '#F196D5', '#A4DDC4', '#589C39', '#18AA59', '#B340C0', '#DFF678', '#5537F4', '#2C866B', '#AB4204',
 '#BB62C5', '#4BA671', '#EAB624', '#3DC26E', '#2748EA', '#328E0B', '#FC046B', '#D18D16', '#04D0EC', '#93AEDE', '#BBDF50', '#23A42B', '#114848']








