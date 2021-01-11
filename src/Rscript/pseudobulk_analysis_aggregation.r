# input
bam_dir <- commandArgs(trailingOnly = T)[1]
pseudoBulk_dir <- commandArgs(trailingOnly = T)[2]
samtoolsbin <- commandArgs(trailingOnly = T)[3]
macs2bin <- commandArgs(trailingOnly = T)[4]
deeptoolsbin <- commandArgs(trailingOnly = T)[5]
genome <- commandArgs(trailingOnly = T)[6]
tabix_path <- commandArgs(trailingOnly = T)[7]
bash_function_dir <- commandArgs(trailingOnly = T)[8]
cores <- commandArgs(trailingOnly = T)[9]

exbin <- commandArgs(trailingOnly = T)[10]
rbin <- commandArgs(trailingOnly = T)[11]
macs2peaktype <- commandArgs(trailingOnly = T)[12]
peak_caller <- commandArgs(trailingOnly = T)[13]

# set the macs2 genome and deeptools effective size
# effectiveGenomeSize https://deeptools.readthedocs.io/en/latest/content/feature/effectiveGenomeSize.html
if(genome=="hg38"){
    macs2_genome <- "hs"
    eGenomeSize <- 2913022398
} else if (genome=="hg19"){
    macs2_genome <- "hs"
    eGenomeSize <- 2864785220
} else if (genome=="mm10"){
    macs2_genome <- "mm"
    eGenomeSize <- 2652783500
} else if (genome=="mm9"){
    macs2_genome <- "mm"
    eGenomeSize <- 2620345972
}

setwd(pseudoBulk_dir)
cat(date(), paste("Processing the group", "[[ groups_aggregation ]]",  "...\n"))
dir.create("groups_aggregation", showWarnings = FALSE)
scbam_dir <- paste0(pseudoBulk_dir, "/", "groups_aggregation/single_cells")
scPS_dir <- paste0(pseudoBulk_dir, "/", "groups_aggregation/pseudo_bulk_data")
# # absolute path
# scbam_dir <- file.path(normalizePath(scbam_dir))
# scPS_dir <- file.path(normalizePath(scPS_dir))
cat(date(), paste("Build the scbam_dir directory", scbam_dir, "...\n"))
cat(date(), paste("Build the scPS_dir directory", scPS_dir, "...\n"))
dir.create(scbam_dir, showWarnings = FALSE)
dir.create(scPS_dir, showWarnings = FALSE)
experi_name <- "groups_aggregation"

cat(date(), paste(length(grep("*bam$", dir(bam_dir))), "Barcode bam files will be copied and merged ...\n"))
myfile <- dir(bam_dir)
# copy the bam files to subgroup dir
file_copy_result <- file.copy(paste(bam_dir, myfile, sep="/"), paste(scbam_dir, myfile, sep="/"), overwrite=T)
message("[info] single-cell track generating")

system(paste("chmod +x", paste0(bash_function_dir, "/qbed.sh")))
system(paste(paste0(bash_function_dir, "/qbed.sh"), scbam_dir, experi_name, tabix_path, scPS_dir))

cat(date(), paste("Pseudobulk bam file", paste0("groups_aggregation_pseudo_sort.bam"), "will be generated"))

if(length(grep("*bam$", dir(scbam_dir), value=T))>1000){
    temp_job_num <- 50
    #message("[info] ", temp_job_num, " temp bam files will be generated")
    mm <- as.integer(seq(1, length(grep("*bam$", dir(scbam_dir), value=T)), length.out=temp_job_num))
    temp_idx_mat <- data.frame(mm[-length(mm)], c((mm-1)[-c(1, length(mm))], mm[length(mm)]))
    for (i in 1:nrow(temp_idx_mat)){
        temp_bam <- paste0(scbam_dir, "/", grep("*bam$", dir(scbam_dir), value=T)[temp_idx_mat[i, 1]:temp_idx_mat[i, 2]], collapse=" ")
        system(paste(paste0(samtoolsbin, "/samtools merge -f"), paste0(scbam_dir, "/", "groups_aggregation", i, "_pseudo.temp"), temp_bam))
        }
    system(paste(paste0(samtoolsbin, "/samtools merge -f"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo.bam"), paste0(scbam_dir, "/*.temp \n")))
    } else {
        message("Merging the individual bam files")
    system(paste(paste0(samtoolsbin, "/samtools merge -f"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo.bam"), paste0(scbam_dir, "/*.bam \n")))
}

system(paste(paste0(samtoolsbin, "/samtools sort -o"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), paste("-@", cores), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo.bam \n")))
system(paste(paste0(samtoolsbin, "/samtools index"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam")))


    # peak calling
    if (peak_caller=="macs2"){
        # macs2
        macs2dir <- paste0(scPS_dir, "/macs2")
        dir.create(macs2dir, showWarnings = FALSE)
        cat(date(), paste("Peak calling for Pseudo-bulk bam file",  paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), "...\n"))
        #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", "groups_aggregation", "macs2.log")))

        if (macs2peaktype=="narrow"){
            message(date(), paste(" [[ MACS2 narrow ]] peaks will be called"))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "_macs2_narrow.log")))
        } else if(macs2peaktype=="broad"){
            message(date(), paste(" [[ MACS2 narrow ]] && [[ MACS2 broad ]] peaks will be called"))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE --broad -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all >", paste0(macs2dir, "/", "groups_aggregation", "_macs2_broad.log")))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "_macs2_narrow.log")))
        }
    } else if (peak_caller=="SEACR"){
        # macs2
        macs2dir <- paste0(scPS_dir, "/macs2")
        dir.create(macs2dir, showWarnings = FALSE)
        cat(date(), paste("Peak calling (MACS2) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), "...\n"))
        #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", "groups_aggregation", "macs2.log")))
        if (macs2peaktype=="narrow"){
            message(date(), paste(" [[ MACS2 narrow ]] peaks will be called"))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "_macs2_narrow.log")))
        } else if(macs2peaktype=="broad"){
            message(date(), paste(" [[ MACS2 narrow ]] && [[ MACS2 broad ]] peaks will be called"))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE --broad -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "_macs2_broad.log")))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "_macs2_narrow.log")))
        }
        # seacr
        seacrdir <- paste0(scPS_dir, "/SEACR")
        dir.create(seacrdir, showWarnings = FALSE)
        message(date(), paste(" [[ SEACR ]] peaks will be called"))
        cat(date(), paste("Peak calling (SEACR) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), "...\n"))
        system(paste("awk '$4!=0'", paste0(macs2dir, "/", "groups_aggregation", "_treat_pileup.bdg"), ">", paste0(seacrdir, "/", "groups_aggregation", "_seacr.bdg")))
        system(paste(paste0(exbin, "/SEACR_1.1.sh"), paste0(seacrdir, "/", "groups_aggregation", "_seacr.bdg"), "0.01 non stringent", paste0(seacrdir, "/", "groups_aggregation"), rbin))
        
    } else {
        message(date(), " Please check your parameter $peak_caller! Will use macs2 narrow peak mode as default.")
        # macs2
        macs2dir <- paste0(scPS_dir, "/macs2")
        dir.create(macs2dir, showWarnings = FALSE)
        cat(date(), paste("Peak calling (macs2) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), "...\n"))
        #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", "groups_aggregation", "macs2.log")))
        message(date(), paste(" ", macs2peaktype, " peaks will be called"))
        system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", "groups_aggregation", "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", "groups_aggregation", "macs2_narrow.log")))
        }
    
    # bam2bw
    cat(date(), paste("Peak calling finished ...\n"))
    cat(date(), paste("Generating the bw file...\n"))
    cat(date(), paste("CPM normalization method will be applied...\n"))
    cat(date(), paste("binsize: 10 \n"))
    cat(date(), paste("Multiple threads will be used: ", cores, " \n"))
    system(paste(paste0(deeptoolsbin, "/bamCoverage"), "--bam", paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bam"), "-o",  paste0(scPS_dir, "/", "groups_aggregation", "_pseudo_sort.bw"), "-p", cores, "--binSize 10 --ignoreDuplicates --normalizeUsing CPM --effectiveGenomeSize", eGenomeSize, paste0("> ", macs2dir, "/", "groups_aggregation_bam2wig.log 2>&1")))




