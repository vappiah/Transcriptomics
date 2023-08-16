# ATTENTION !!!!!!!!!!!!!!!!!!!!!

This command '''keep <- rowSums(counts(dds)>10)'''

Shoud be written as '''keep <- rowSums(counts(dds)>10) >= min(table(sample_info$Treatment)) '''
