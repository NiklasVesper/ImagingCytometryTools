def phenotype_row(df):
    cell_types_and_states = []  # empty list for all cell types

    for index, cell in df.iterrows():
        cell_type = []

        CD45 = []
        if cell['Intensity_MeanIntensity_CD45_Cell'] >= 0.15:
            CD45.append('CD45+')

            CD4 = []
            if cell['Intensity_MeanIntensity_CD4_Cell'] > 0.15:
                CD4.append('CD4+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD4.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD4.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    CD4.append('CXCR5+')
                if cell['Intensity_MeanIntensity_CXCR6_Cell'] > 0.5:
                    CD4.append('CXCR6+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    CD4.append('Tbet+')
                if cell['Intensity_MeanIntensity_TCF1_Nucleus'] > 0.5:
                    CD4.append('TCF1+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    CD4.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD4.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD4.append('CD39+')
                if cell['Intensity_MeanIntensity_CD45RO_Cell'] > 0.4:
                    CD4.append('CD45RO+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD4.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD4.append('FAK+')

            FoxP3 = []

            if cell['Intensity_MeanIntensity_FoxP3_Nucleus'] > 0.2 and cell['Intensity_MeanIntensity_CD4_Cell'] > 0.15:
                FoxP3.append('FoxP3+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    FoxP3.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    FoxP3.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    FoxP3.append('CXCR5+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    FoxP3.append('Tbet+')
                if cell['Intensity_MeanIntensity_TCF1_Nucleus'] > 0.5:
                    FoxP3.append('TCF1+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    FoxP3.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    FoxP3.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    FoxP3.append('CD39+')
                if cell['Intensity_MeanIntensity_CD45RO_Cell'] > 0.4:
                    FoxP3.append('CD45RO+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    FoxP3.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    FoxP3.append('FAK+')

            if len(FoxP3) >= 1:
                cell_type.append(FoxP3)
            elif len(CD4) >= 1:
                cell_type.append(CD4)

            CD8 = []
            if cell['Intensity_MeanIntensity_CD8_Cell'] > 0.5:
                CD8.append('CD8+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD8.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD8.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    CD8.append('CXCR5+')
                if cell['Intensity_MeanIntensity_CXCR6_Cell'] > 0.5:
                    CD8.append('CXCR6+')
                if cell['Intensity_MeanIntensity_GranzymeB_Cell'] > 0.5:
                    CD8.append('GranzymeB+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    CD8.append('Tbet+')
                if cell['Intensity_MeanIntensity_TCF1_Nucleus'] > 0.5:
                    CD8.append('TCF1+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    CD8.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD8.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD8.append('CD39+')
                if cell['Intensity_MeanIntensity_CD45RO_Cell'] > 0.4:
                    CD8.append('CD45RO+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD8.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD8.append('FAK+')

            if len(CD8) >= 1:
                cell_type.append(CD8)

            CD20 = []
            if cell['Intensity_MeanIntensity_CD20_Cell'] > 0.35:
                CD20.append('CD20+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD20.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD20.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    CD20.append('CXCR5+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    CD20.append('Tbet+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    CD20.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD20.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD20.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD20.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD20.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD20.append('FAK+')

            if len(CD20) >= 1:
                cell_type.append(CD20)

            CD15 = []
            if cell['Intensity_MeanIntensity_CD15_Cell'] > 0.4:
                CD15.append('CD15+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD15.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD15.append('PD1+')
                if cell['Intensity_MeanIntensity_GranzymeB_Cell'] > 0.5:
                    CD15.append('GranzymeB+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD15.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD15.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD15.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD15.append('FAK+')


            if len(CD15) >= 1:
                cell_type.append(CD15)

            CD68 = []
            if cell['Intensity_MeanIntensity_CD68_Cell'] > 0.6:
                CD68.append('CD68+')

                if cell['Intensity_MeanIntensity_CD204_Cell'] > 0.4:
                    CD68.append('CD204+')
                if cell['Intensity_MeanIntensity_CD163_Cell'] > 0.3:
                    CD68.append('CD163+')
                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD68.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD68.append('PD1+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD68.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD68.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD68.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD68.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD68.append('FAK+')

            if len(CD68) >= 1:
                cell_type.append(CD68)

            CD11c = []
            if cell['Intensity_MeanIntensity_CD11c_Cell'] > 0.4:
                CD11c.append('CD11c+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD11c.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD11c.append('PD1+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 1:
                    CD11c.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD11c.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD11c.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD11c.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD11c.append('FAK+')


            if len(CD11c) >= 1:
                cell_type.append(CD11c)

        elif cell['Intensity_MeanIntensity_CD45_Cell'] < 0.15:

            CD45.append('CD45-')

            FoxP3 = []
            if cell['Intensity_MeanIntensity_FoxP3_Nucleus'] > 0.3 and cell['Intensity_MeanIntensity_CD4_Cell'] > 0.15:
                FoxP3.append('FoxP3+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    FoxP3.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    FoxP3.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    FoxP3.append('CXCR5+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    FoxP3.append('Tbet+')
                if cell['Intensity_MeanIntensity_TCF1_Nucleus'] > 0.5:
                    FoxP3.append('TCF1+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    FoxP3.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    FoxP3.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    FoxP3.append('CD39+')
                if cell['Intensity_MeanIntensity_CD45RO_Cell'] > 0.4:
                    FoxP3.append('CD45RO+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    FoxP3.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    FoxP3.append('FAK+')

            if len(FoxP3) >= 1:
                cell_type.append(FoxP3)
                CD45.append('CD45+')

            CD8 = []
            if cell['Intensity_MeanIntensity_CD8_Cell'] > 0.6:
                CD8.append('CD8+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD8.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD8.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    CD8.append('CXCR5+')
                if cell['Intensity_MeanIntensity_CXCR6_Cell'] > 0.5:
                    CD8.append('CXCR6+')
                if cell['Intensity_MeanIntensity_GranzymeB_Cell'] > 0.5:
                    CD8.append('GranzymeB+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    CD8.append('Tbet+')
                if cell['Intensity_MeanIntensity_TCF1_Nucleus'] > 0.5:
                    CD8.append('TCF1+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    CD8.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD8.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD8.append('CD39+')
                if cell['Intensity_MeanIntensity_CD45RO_Cell'] > 0.4:
                    CD8.append('CD45RO+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD8.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD8.append('FAK+')

            if len(CD8) >= 1:
                cell_type.append(CD8)
                CD45.append('CD45+')

            CD20 = []
            if cell['Intensity_MeanIntensity_CD20_Cell'] > 0.4:
                CD20.append('CD20+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD20.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD20.append('PD1+')
                if cell['Intensity_MeanIntensity_CXCR5_Cell'] > 1:
                    CD20.append('CXCR5+')
                if cell['Intensity_MeanIntensity_Tbet_Nucleus'] > 0.2:
                    CD20.append('Tbet+')
                if cell['Intensity_MeanIntensity_TOX_Nucleus'] > 0.25:
                    CD20.append('TOX+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD20.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD20.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD20.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD20.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD20.append('FAK+')

            if len(CD20) >= 1:
                cell_type.append(CD20)
                CD45.append('CD45+')

            CD15 = []
            if cell['Intensity_MeanIntensity_CD15_Cell'] > 0.45:
                CD15.append('CD15+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD15.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD15.append('PD1+')
                if cell['Intensity_MeanIntensity_GranzymeB_Cell'] > 0.5:
                    CD15.append('GranzymeB+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD15.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD15.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD15.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD15.append('FAK+')

            if len(CD15) >= 1:
                cell_type.append(CD15)
                CD45.append('CD45+')

            CD68 = []
            if cell['Intensity_MeanIntensity_CD68_Cell'] > 0.6:
                CD68.append('CD68+')

                if cell['Intensity_MeanIntensity_CD204_Cell'] > 0.4:
                    CD68.append('CD204+')
                if cell['Intensity_MeanIntensity_CD163_Cell'] > 0.3:
                    CD68.append('CD163+')
                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD68.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD68.append('PD1+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 0.5:
                    CD68.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD68.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD68.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD68.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD68.append('FAK+')

            if len(CD68) >= 1:
                cell_type.append(CD68)
                CD45.append('CD45+')

            CD11c = []
            if cell['Intensity_MeanIntensity_CD11c_Cell'] > 0.45:
                CD11c.append('CD11c+')

                if cell['Intensity_MeanIntensity_CD38_Cell'] > 0.25:
                    CD11c.append('CD38+')
                if cell['Intensity_MeanIntensity_PD1_Cell'] > 0.2:
                    CD11c.append('PD1+')
                if cell['Intensity_MeanIntensity_Tim3_Cell'] > 1:
                    CD11c.append('Tim3+')
                if cell['Intensity_MeanIntensity_CD39_Cell'] > 0.25:
                    CD11c.append('CD39+')
                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.3:
                    CD11c.append('ATF6+')
                if cell['Intensity_MeanIntensity_HLADR_Cell'] > 0.5:
                    CD11c.append('HLADR+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 1.5:
                    CD11c.append('FAK+')

            if len(CD11c) >= 1:
                cell_type.append(CD11c)
                CD45.append('CD45+')

            Hepatocyte = []
            if cell['Intensity_MeanIntensity_Ecadherin_Cell'] > 0.5 and cell['Intensity_MeanIntensity_Bcatenin_Cell'] > 0.5 and cell['Intensity_MeanIntensity_FoxP3_Nucleus'] <= 0.3 and cell['Intensity_MeanIntensity_CD4_Cell'] <= 0.15 and cell['Intensity_MeanIntensity_CD8_Cell'] <= 0.6 and cell['Intensity_MeanIntensity_CD20_Cell'] <= 0.4 and cell['Intensity_MeanIntensity_CD15_Cell'] <= 0.45 and cell['Intensity_MeanIntensity_CD68_Cell'] <= 0.6 and cell['Intensity_MeanIntensity_CD11c_Cell'] <= 0.45:
                Hepatocyte.append('ECadherin+ BCatenin+')

                if cell['Intensity_MeanIntensity_ATF6_Nucleus'] > 0.2:
                    Hepatocyte.append('ATF6+')
                if cell['Intensity_MeanIntensity_FAK_Nucleus'] > 2:
                    Hepatocyte.append('FAK+')

            if len(Hepatocyte) >= 1:
                cell_type.append(Hepatocyte)

        CD45_final = list(set(CD45))

        if len(CD45_final) == 1:
            cell_type.append(CD45_final)

        if len(CD45_final) > 1:
            cell_type.append(['CD45+'])

        cell_types_and_states.append(cell_type)

    df['cell_types_and_states'] = cell_types_and_states




meta_data_list_1 = ['example_name_1_subcellular_and_neighborhood_5',
                    'example_name_2_subcellular_and_neighborhood_5',
                    'example_name_3_subcellular_and_neighborhood_5',
                    'example_name_4_subcellular_and_neighborhood_5',
                    'example_name_5_subcellular_and_neighborhood_5']

meta_data_list_2 = ['example_name_6_subcellular_and_neighborhood_5',
                  'example_name_7_subcellular_and_neighborhood_5',
                  'example_name_8_subcellular_and_neighborhood_5',
                  'example_name_9_subcellular_and_neighborhood_5',
                  'example_name_10_subcellular_and_neighborhood_5']

sample_parings={'example_name_1_subcellular_and_neighborhood_5':'example_name_6_subcellular_and_neighborhood_5',
                'example_name_2_subcellular_and_neighborhood_5':'example_name_7_subcellular_and_neighborhood_5',
                'example_name_3_subcellular_and_neighborhood_5':'example_name_8_subcellular_and_neighborhood_5',
                'example_name_4_subcellular_and_neighborhood_5':'example_name_9_subcellular_and_neighborhood_5',
                'example_name_5_subcellular_and_neighborhood_5':'example_name_10_subcellular_and_neighborhood_5'}

