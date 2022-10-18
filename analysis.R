##
library(NetBID2)
source('functions_haplotypeView.R')
##
WGS_block <- read_block('data/phasedWGS.plink.blocks.blocks.det','data/phasedWGS.plink.map')
carrier_block <- read_block('data/phased.impute2.blocks.det',info_file='data/phased.impute2.info')
## basic process
carrier <- read_ped('data/phased.impute2.ped',info_file='data/phased.impute2.info')
CES <- read_ped('data/phasedCES.impute2.ped',info_file='data/phasedCES.impute2.info')
eset_carrier <- merge_ped(list(carrier=carrier)) ## 815*30
eset_carrier <- combine_pos_ped(eset_carrier) ## 216*30
eset_both <- merge_ped(list(CES=CES,carrier=carrier))
eset_both <- combine_pos_ped(eset_both) ## 543*3918
eset_carrier <- filter_sample(eset_carrier,use_sample=setdiff(pData(eset_carrier)$ori_sample_name,c('R17001732LU01','P19016041LU01')))
eset_both <- filter_sample(eset_both,use_sample=setdiff(pData(eset_both)$ori_sample_name,c('R17001732LU01','P19016041LU01')))
##
eset11 <- filter_pos(eset_carrier,thre=0.001)
eset12 <- filter_pos(eset_carrier,thre=0.001,use_pos=unlist(lapply(WGS_block,function(x)x$mark_pos)),remove_identical=T)
eset_both11 <- filter_pos(eset_both,thre=0.001)
eset_both12 <- filter_pos(eset_both,thre=0.001,use_pos=unlist(lapply(WGS_block,function(x)x$mark_pos)))
eset_random <- random_select_eset(eset_both12,n=100,use_sample=rownames(pData(eset11)))
## visualization mark_pos=c(51375699,51379280)
draw_pos2sample_pdf('t1.pdf',eset12,remove_identical=TRUE,mark_pos=c(51375699,51379280),
                    order_strategy='name',mark_block=WGS_block,point_strategy='text',
                    top_hap_number=6,only_mark_tag_SNP=TRUE)
draw_pos2sample_pdf('t2.pdf',eset_random,remove_identical=F,mark_pos=c(51375699,51379280),
                    order_strategy='name',mark_block=WGS_block,point_strategy='text',
                    top_hap_number=6,only_mark_tag_SNP=TRUE)
## hap assoc test ?
eset_both12 <- filter_pos(eset_both12,thre=0.001)
eset_hap <- get_block_hap(eset_both12,mark_block=WGS_block,remove_miss=T)
res <- test_hap_assoc(eset_hap,group='group')
res1 <- test_hap_assoc(eset_both12,group='group')



