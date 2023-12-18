import pandas as pd


'''
    When taking data from TCGA, select columns (three horizontal bars): DNA Change, #Affected Case in Cohort, and Impact
    Select top 100 rows
'''
colon_data = pd.read_csv('cancer_data/colon.tsv', sep='\t')
rectum_data = pd.read_csv('cancer_data/rectum.tsv', sep='\t')
location1 = 'Colon'
location2 = 'Rectum'


""" Helper Methods """

def to_float(fraction: str):
    if isinstance(fraction, str):
        a, b = fraction.replace(',', '').split('/')
        return round(float(a) / float(b), 4)
    else:
        return fraction


'''------------------------------------------------------------------------------------'''

'''
    Join data from two dataframes into one easy-to-read dataframe
    Parameters: 2 dataframes
'''
def merge_data(df_1, df_2):
    combined_data = pd.merge(df_1, df_2, on='DNA Change', how='outer')
    combined_data.drop(['Impact_x'], axis=1, inplace=True)  
    combined_data.rename(columns={'# Affected Cases in Cohort_x': f'{location1} Frequency', '# Affected Cases in Cohort_y': f'{location2} Frequency', 'Impact_y': 'Impact', '# Affected Cases Across the GDC_x' : 'Total'}, inplace=True)
    combined_data[['Fraction', f'{location1} Cases']] = combined_data[f'{location1} Frequency'].str.split(',', expand=True)
    combined_data[['Fraction2', f'{location2} Cases']] = combined_data[f'{location2} Frequency'].str.split(',', expand=True)
    combined_data['Total'] = combined_data['Total'].combine_first(combined_data['# Affected Cases Across the GDC_y'])
    combined_data['Total Cases'] = combined_data['Total'].apply(lambda x: to_float(x))
    combined_data['Colon Cases'] = combined_data['Colon Cases'].str.rstrip('%').astype('float') / 100.0
    combined_data['Rectum Cases'] = combined_data['Rectum Cases'].str.rstrip('%').astype('float') / 100.0
    combined_data.drop([f'{location1} Frequency', f'{location2} Frequency', 'Fraction', 'Fraction2', '# Affected Cases Across the GDC_y', 'Impact'], axis=1, inplace=True)
    return combined_data

'''------------------------------------------------------------------------------------'''

'''
    Finds mutations that are in one form of cancer but not the other
'''
def find_differences(df):
    data1_mask = df[f'{location1} Cases'].notna() & df[f'{location2} Cases'].isna()
    data2_mask = df[f'{location2} Cases'].notna() & df[f'{location1} Cases'].isna()

    data1_mutations = df[data1_mask]
    data2_mutations = df[data2_mask]

    return ([data1_mutations['DNA Change'], data2_mutations['DNA Change']])

'''------------------------------------------------------------------------------------'''

'''
    Find number of mutations in only one of the cancer forms
'''
def find_num_differences(df):
    return f"{location1}: {len(find_differences(df)[0])}\n{location2}: {len(find_differences(df)[1])}"

'''------------------------------------------------------------------------------------'''

'''
    Find common mutations in two cancer forms
'''
def find_commonalities(df_1, df_2):
    joined_data = pd.merge(df_1, df_2, on='DNA Change', how='inner')
    joined_data.rename(columns={'# Affected Cases in Cohort_x': f'{location1} Frequency', '# Affected Cases in Cohort_y': f'{location2} Frequency', 'Impact_y': 'Impact'}, inplace=True)
    return joined_data

'''------------------------------------------------------------------------------------'''

'''
    Find number of common mutations in two cancer forms
'''
def find_num_commonalities(df_1, df_2):
    return len(find_commonalities(df_1, df_2))

'''------------------------------------------------------------------------------------'''

'''
    Parameters: merged dataframe, DNA Change value
    Finds if a DNA change is found in both cancer forms
    If it finds a match in both, it returns the proportion of both relative to entire cohort
'''

def in_both_sets(df, dna_change):
    index = df.index[df['DNA Change'] == dna_change]
    if not index.empty:
        index = index[0]  # Assuming you want the first matching index if there are multiple
        colon_cases = df.at[index, 'Colon Cases']
        rectum_cases = df.at[index, 'Rectum Cases']
        return colon_cases, rectum_cases
    return False

#Input merged data, cancer location to be compared (capitalize first letter)
def compare_to_total(merged_data, location):
    case = pd.DataFrame()
    case['DNA Change'] = merged_data['DNA Change']
    case[f'{location} to Total'] = merged_data[f'{location} Cases'] - merged_data['Total Cases']
    return case.dropna(subset = [f'{location} to Total']).sort_values(by = ['Colon to Total'], ascending=False)

#Merged Colon and Rectum Data
merged_data = merge_data(colon_data, rectum_data)

print(find_num_differences(merged_data))

print(in_both_sets(merged_data, "chr7:g.140753336A>T"))
