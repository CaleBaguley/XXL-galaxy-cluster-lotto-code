"""
Converts csv list of positions to text file for passing to HSC data request.
"""


import pandas as pd

source_list = 'Spurious sources removed by preprocessing/Sources_removed_by_cleaning.csv'

max_source_count = 100


# Column keys for file
Name_col = 'Index'
Ra_col   = 'EXT_RA'
dec_col  = 'EXT_DEC'

# Size of cutout (degrees)
size = 5

# output files. Set to None to skip
wide_output = 'Spurious sources removed by preprocessing/Sources_removed_by_cleaning_selected.txt'
dud_output = None

# Load data
data = pd.read_csv(source_list)



if(max_source_count is not None and max_source_count < len(data)):
    print("randomly reducing data from {} to {}".format(len(data), max_source_count))
    data = data.sample(n = max_source_count)

print(data)

# Sort by Ra to speed up download
data.sort_values(Ra_col, inplace=True)

# Wide survey request
if(isinstance(wide_output,str)):
    # open file to write to
    with open(wide_output, 'w') as f:
        
        # write first line
        f.write('#?  name  ra  dec  sw  sh  filter  rerun\n')
        
        for i, entry in data.iterrows():
            f.write('  {}_G  {}  {}  {}amin  {}amin  HSC-G  pdr3_wide\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_R  {}  {}  {}amin  {}amin  HSC-R  pdr3_wide\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_I  {}  {}  {}amin  {}amin  HSC-I  pdr3_wide\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_Z  {}  {}  {}amin  {}amin  HSC-Z  pdr3_wide\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))  


# Dud survey request
if(isinstance(dud_output,str)):
    # open file to write to
    with open(dud_output, 'w') as f:
        
        # write first line
        f.write('#?  name  ra  dec  sw  sh  filter  rerun\n')
        
        for i, entry in data.iterrows():
            f.write('  {}_G  {}  {}  {}amin  {}amin  HSC-G  pdr3_dud\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_R  {}  {}  {}amin  {}amin  HSC-R  pdr3_dud\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_I  {}  {}  {}amin  {}amin  HSC-I  pdr3_dud\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
        
        for i, entry in data.iterrows():
            f.write('  {}_Z  {}  {}  {}amin  {}amin  HSC-Z  pdr3_dud\n'.format(str(entry[Name_col]), entry[Ra_col], entry[dec_col], size/2, size/2))
            
            
            
            
            
