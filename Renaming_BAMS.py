import os
import os.path



SAMPLE_TO_SPECIES={}

if os.path.exists('METADATA'):
    META_DATA_FILE=open('METADATA','r')


if os.path.exists('METADATA'):
    META_DATA_FILE.readline()


    #### Read METADATA file, create Dictionary
    

    
    for LINE in META_DATA_FILE:
    

        LINE=LINE.strip().split('\t')
        SAMPLE_NAME=LINE[1]
        SPECIES='_'.join(LINE[5].split(' '))
        SAMPLE_TO_SPECIES[SAMPLE_NAME]=SPECIES



print(SAMPLE_TO_SPECIES)
print('\n')




for V,J in SAMPLE_TO_SPECIES.items():
    
    os.system(f'cp {V}.bam FINISHED/{V}_{J}.bam')
    os.system(f'cp {V}.bam.bai FINISHED/{V}_{J}.bam.bai')