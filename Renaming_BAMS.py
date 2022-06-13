import os
import os.path



BAM_TO_NAME={}

if os.path.exists('METADATA'):
    META_DATA_FILE=open('METADATA','r')


if os.path.exists('METADATA'):
    META_DATA_FILE.readline()


    #### Read METADATA file, create Dictionary
    
    SPECIES_COUNT={}

    
    for LINE in META_DATA_FILE:
    

        LINE=LINE.strip().split('\t')
        SAMPLE=LINE[1]
        SPECIES='_'.join(LINE[5].split(' '))
        NAME=SPECIES
        
        
        if SAMPLE not in BAM_TO_NAME.keys():
            if NAME not in SPECIES_COUNT.keys():
                SPECIES_COUNT[NAME]=1
            else:
                SPECIES_COUNT[NAME]+=1
        
            print(NAME,SAMPLE)
                

        BAM_TO_NAME[SAMPLE]=NAME


print(BAM_TO_NAME)
print('\n')
print(SPECIES_COUNT)


for J in SPECIES_COUNT.keys():
    counter=1
    for K in BAM_TO_NAME.keys():
        
        if BAM_TO_NAME[K]==J:
            BAM_TO_NAME[K]=BAM_TO_NAME[K]+'_'+str(counter)
            counter+=1

print(BAM_TO_NAME)

for V in BAM_TO_NAME:
    
    os.system(f'mv FINISHED/{V}.bam FINISHED/{V}_{BAM_TO_NAME[V]}.bam')
    os.system(f'mv FINISHED/{V}.bam.bai FINISHED/{V}_{BAM_TO_NAME[V]}.bam.bai')