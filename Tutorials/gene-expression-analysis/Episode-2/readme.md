# ATTENTION !!!!!!!!!!!!!!!!!!!!!
# THIS SHOULD BE CORRECTED WHEN WATCHING THE TUTORIAL !!!!!!!!!!!!!!!!!!!!!!!11

This command '''keep <- rowSums(counts(dds)>10)'''

Shoud be written as '''keep <- rowSums(counts(dds)>10) >= min(table(sample_info$Treatment)) '''
