library(ShortRead)
library(ggplot2)
library(dplyr)
library(Rsubread)
library(cowplot)
library(ggplot2)
library(gridExtra)
# Set paths

setwd("/home/data2/Swapna/RNAQuantSeq")
fastq_dirs <- "/home/data2/Swapna/RNAQuantSeq/fastQfiles"
sinkpath  <- "/home/data2/Swapna/RNAQuantSeq/bam_files"
ref_fa <- "/home/data2/Swapna/RNAQuantSeq/ref_dir/Homo_sapiens.GRCh38.dna_sm.chromosome.1.fa"
ref_dir <- "/home/data2/Swapna/RNAQuantSeq/ref_dir/"
gtf.file <- file.path(ref_dir,"Homo_sapiens.GRCh38.84.gtf")
sqlite_file <- 'Homo_sapiens.GRCh38.84.sqlite'
sqlite_path <- file.path(ref_dir, sqlite_file)
workflow_stats <- list()

# Quality check plots fucntion

qa_plots <- function(shortfq_obj) {
    x <- ShortRead::qa(shortfq_obj, lane='')
    nuc_plot <- ggplot(x[['perCycle']]$baseCall, aes(x=Cycle, y=Count, colour=Base)) + geom_line() + cowplot::theme_cowplot()
    qual_plot <- x[['perCycle']]$quality %>% group_by(Cycle) %>% summarise(QualScore=sum(Score*Count)/sum(Count)) %>% ungroup() %>% 
        ggplot(aes(x=Cycle, y=QualScore)) + geom_line() + cowplot::theme_cowplot()
    size_plot <- ggplot(data.frame(frag_size=width(shortfq_obj)), aes(x=frag_size)) + 
        geom_bar(fill='pink', colour='black') + cowplot::theme_cowplot()
    base_plot <- ggplot(data.frame(idx=1:50, dat=as.character(sread(shortfq_obj)[1:50]))) + 
        geom_text(aes(x=1, y=idx, label=dat), size=rel(2.5), family='mono', hjust=0) + xlim(1,1.8) +
        theme_bw() + theme(axis.text.x=element_text(size=0)) + xlab('')
    return(list(nuc_plot=nuc_plot, qual_plot=qual_plot, size_plot=size_plot, base_plot=base_plot))
}

# Import data

fq_fn <- list.files(fastq_dirs, 'fastq.gz')
fq_path <- file.path(fastq_dirs, fq_fn)
fp <- dput(as.character((fq_path)))

fq_data_list <- list()
for(i in 1:3){
fq_data_list[[i]] <- yield(FastqSampler(fp[[i]]))
}

print( fq_data_list)

# QC plots 
## Generates one .pdf per sample.

qp <- list()
for ( i in 1:3){
ggsave(marrangeGrob(qa_plots(fq_data_list[[i]]), nrow = 3, ncol =1), file=paste0(i,"_1.pdf"))
}
 
# Remove the 5’ Adapter (first 12 bps)

fq <- fq_data_list

fq_list <- lapply(fq,ShortRead::narrow,start = 13)
names(fq_list) <- gsub(".fastq.gz","",fq_fn[1:3])

# confirm that the adapters have been removed

qp <- list()
for ( i in 1:3){
ggsave(marrangeGrob(qa_plots(fq_list[[i]]), nrow = 3, ncol =1), file=paste0(i,"_2.pdf"))
}

# Remove poor quality data at the 3’ end of the read

fq_list <- lapply(fq_list, ShortRead::trimTailw, k = 6, a="4", halfwidth=6)
for( i in 1:3){
fq_list[[i]] <- fq_list[[i]][width(fq_list[[i]]) >= 36]
print(length(fq_list[[i]]))
}

# plot

for ( i in 1:3){
ggsave(marrangeGrob(qa_plots(fq_list[[i]]), nrow = 3, ncol =1), file=paste0(i,"_3.pdf"))
}
 
# Trim polyA tail


narrow_ranges_list <- lapply(lapply(fq_list,sread), ShortRead::trimEnds, right=TRUE, "A", relation="==", ranges=TRUE)

fq_list1 <- list()
for( i in 1:3){
fq_list1[[i]] <- ShortRead::narrow(fq_list[[i]], start(narrow_ranges_list[[i]]), end(narrow_ranges_list[[i]]))
}

fq_list2 <- list()
for( i in 1:3){
fq_list2[[i]] <- fq_list1[[i]][width(fq_list1[[i]]) >= 36]
print(length(fq_list2[[i]]))
}

# plot

for ( i in 1:3){
ggsave(marrangeGrob(qa_plots(fq_list2[[i]]), nrow = 3, ncol =1), file=paste0(i,"_4.pdf"))
}
 
# Trim embedded polyA

polyApos <- list()
polyAclip_idx <- list()
polyAclip <- list()
for ( i in 1:3) {
polyApos[[i]] <- regexpr(pattern= 'AAAAAAAAAA', sread(fq_list2[[i]], fixed=TRUE) %>% unlist())
print(polyApos[[i]])
polyAclip_idx[[i]] <- which(polyApos[[i]] >= 0)

polyAclip[[i]] <- width(fq_list2[[i]])
polyAclip[[i]][polyAclip_idx[[i]]] <- polyApos[[i]][polyAclip_idx[[i]]] 
fq_list2[[i]] <- narrow(fq_list2[[i]], end=polyAclip[[i]])
fq_list2[[i]] <- fq_list2[[i]] [width(fq_list2[[i]] ) >= 36]
length(fq_list2)[[i]]
}

## returned -1 for all the 3; no embedded poly A found.


# create reference

# buildindex(basename="reference_index",reference=ref_fa)

# alignment 

fastQsub_path <- "/home/data2/Swapna/RNAQuantSeq/trimmedFastq_files/"
fq_tr <- list.files(fastQsub_path, 'trimmedFastq')

bam_path <- "/home/data2/Swapna/RNAQuantSeq/bam_files"

for( i in 1:3){
fq_tr_path[i] <- file.path(fastQsub_path, fq_tr[i])
writeFastq(fq_list[[i]], paste(as.character(fastQsub_path),names(fq_list)[i],sep="_"))


align (readfile1=fq_tr_path[i], index=file.path(ref_dir,"reference_index"),
      input_format='gzFASTQ',
      output_format = "BAM",
      output_file=paste(as.character(fq_tr_path)[i],"subread_BAM",  sep= "."),   
      nthreads=8)

}


