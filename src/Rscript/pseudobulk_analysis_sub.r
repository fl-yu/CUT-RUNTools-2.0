# input
cell_anno <- commandArgs(trailingOnly = T)[1]
bam_dir <- commandArgs(trailingOnly = T)[2]
pseudoBulk_dir <- commandArgs(trailingOnly = T)[3]
samtoolsbin <- commandArgs(trailingOnly = T)[4]
macs2bin <- commandArgs(trailingOnly = T)[5]
deeptoolsbin <- commandArgs(trailingOnly = T)[6]
genome <- commandArgs(trailingOnly = T)[7]
tabix_path <- commandArgs(trailingOnly = T)[8]
bash_function_dir <- commandArgs(trailingOnly = T)[9]
cores <- commandArgs(trailingOnly = T)[10]
exbin <- commandArgs(trailingOnly = T)[11]
rbin <- commandArgs(trailingOnly = T)[12]
macs2peaktype <- commandArgs(trailingOnly = T)[13]
peak_caller <- commandArgs(trailingOnly = T)[14]

# read new meta_table; no rownames
anno_meta <- read.table(cell_anno, sep="\t", header=T)

target_group <- unique(anno_meta[, "leiden_cluster"])
group_num <- length(target_group)

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
output_all_group <- TRUE
if (group_num==1) {
    output_separated_group <- FALSE
    stop("[error] Make sure number of subgroup is more than 1")
} else {
    cat(date(), paste("[info]", group_num, "groups were included for the analysis: \n"))
    for (i in 1:group_num){

        cat(date(), paste("Processing the group", "[[", target_group[i], "]]",  "...\n"))
        dir.create(paste0(paste0("group_", target_group[i])), showWarnings = FALSE)
        scbam_dir <- paste0(paste0(pseudoBulk_dir, "/", paste0("group_", target_group[i])), "/single_cells")
        scPS_dir <- paste0(paste0(pseudoBulk_dir, "/", paste0("group_", target_group[i])), "/pseudo_bulk_data")

        # # absolute path > will change to relative path when use system(), don't know why
        # scbam_dir <- file.path(normalizePath(scbam_dir))
        # scPS_dir <- file.path(normalizePath(scPS_dir))

        cat(date(), paste("Build the scbam_dir directory", scbam_dir, "...\n"))
        cat(date(), paste("Build the scPS_dir directory", scPS_dir, "...\n"))
        dir.create(scbam_dir, showWarnings = FALSE)
        dir.create(scPS_dir, showWarnings = FALSE)
        
        idx <- anno_meta[, "leiden_cluster"] == target_group[i]
        cat(date(), paste(sum(idx), "files will be copied and merged ...\n"))
        single_bam_file <- paste0(bam_dir, "/", anno_meta$cell_name[idx])
        single_bed_file <- gsub(".bam",  ".bed", single_bam_file)
        # copy the bam files to subgroup dir
        file.copy(single_bam_file, scbam_dir)
        file.copy(single_bed_file, scbam_dir)
        message("[info] single-cell track generating")
        system(paste("chmod +x", paste0(bash_function_dir, "/qbed.sh")))
        experi_name <- paste0(paste0("group_", target_group[i]))
        system(paste(paste0(bash_function_dir, "/qbed.sh"), scbam_dir, experi_name, tabix_path, scPS_dir))
        # merge/sort/index/bigwig
        cat(date(), paste("Pseudobulk bam file", paste0(paste0(paste0("group_", target_group[i])), "_pseudo.bam"), "will be generated && sort && indexed in the scPS_dir directory", scPS_dir, "...\n"))
        # samtools sort -o sorted.bam initial.bam
        # samtools index aln.sorted.bam
        system(paste(paste0(samtoolsbin, "/samtools merge"), paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo.bam"), paste0(scbam_dir, "/*.bam \n")))
        system(paste(paste0(samtoolsbin, "/samtools sort -o"), paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo_sort.bam"), paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo.bam \n")))
        system(paste(paste0(samtoolsbin, "/samtools index"), paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo_sort.bam")))
        
        # peak calling
        if (peak_caller=="macs2"){
            # macs2
            macs2dir <- paste0(scPS_dir, "/macs2")
            dir.create(macs2dir, showWarnings = FALSE)
            cat(date(), paste("Peak calling for Pseudo-bulk bam file",  paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"), "...\n"))
            #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "macs2.log")))

            if (macs2peaktype=="narrow"){
                message(date(), paste(" [[ MACS2 narrow ]] peaks will be called"))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_narrow.log")))
            } else if(macs2peaktype=="broad"){
                message(date(), paste(" [[ MACS2 narrow ]] && [[ MACS2 broad ]] peaks will be called"))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE --broad -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all >", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_broad.log")))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_narrow.log")))
            }
        } else if (peak_caller=="SEACR"){
            # macs2
            macs2dir <- paste0(scPS_dir, "/macs2")
            dir.create(macs2dir, showWarnings = FALSE)
            cat(date(), paste("Peak calling (MACS2) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"), "...\n"))
            #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "macs2.log")))
            if (macs2peaktype=="narrow"){
                message(date(), paste(" [[ MACS2 narrow ]] peaks will be called"))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_narrow.log")))
            } else if(macs2peaktype=="broad"){
                message(date(), paste(" [[ MACS2 narrow ]] && [[ MACS2 broad ]] peaks will be called"))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE --broad -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_broad.log")))
                system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_macs2_narrow.log")))
            }
            # seacr
            seacrdir <- paste0(scPS_dir, "/SEACR")
            dir.create(seacrdir, showWarnings = FALSE)
            message(date(), paste(" [[ SEACR ]] peaks will be called"))
            cat(date(), paste("Peak calling (SEACR) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"), "...\n"))
            system(paste("awk '$4!=0'", paste0(macs2dir, "/", paste0("group_", target_group[i]), "_treat_pileup.bdg"), ">", paste0(seacrdir, "/", paste0("group_", target_group[i]), "_seacr.bdg")))
            system(paste(paste0(exbin, "/SEACR_1.1.sh"), paste0(seacrdir, "/", paste0("group_", target_group[i]), "_seacr.bdg"), "0.01 non stringent", paste0(seacrdir, "/", paste0("group_", target_group[i])), rbin))
            
        } else {
            message(date(), " Please check your parameter $peak_caller! Will use macs2 narrow peak mode as default.")
            # macs2
            macs2dir <- paste0(scPS_dir, "/macs2")
            dir.create(macs2dir, showWarnings = FALSE)
            cat(date(), paste("Peak calling (macs2) for Pseudo-bulk bam file",  paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"), "...\n"))
            #system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup 1 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "macs2.log")))
            message(date(), paste(" ", macs2peaktype, " peaks will be called"))
            system(paste(paste0(macs2bin, "/macs2 callpeak -t"), paste0(scPS_dir, "/", paste0("group_", target_group[i]), "_pseudo_sort.bam"),  "-g", macs2_genome, "-f BAMPE -n", paste0("group_", target_group[i]), "--outdir", macs2dir, "-q 0.01 -B --SPMR --keep-dup all 2>", paste0(macs2dir, "/", paste0("group_", target_group[i]), "macs2_narrow.log")))
            }
        cat(date(), paste("Peak calling finished ...\n"))
        cat(date(), paste("Generating the bw file...\n"))
        cat(date(), paste("CPM normalization method will be applied...\n"))
        cat(date(), paste("binsize: 10 \n"))
        cat(date(), paste("Multiple threads will be used: ", cores, " \n"))
        # bam2bw   
        system(paste(paste0(deeptoolsbin, "/bamCoverage"), "--bam", paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo_sort.bam"), "-o",  paste0(scPS_dir, "/", paste0(paste0("group_", target_group[i])), "_pseudo_sort.bw"), "-p", cores, "--binSize 10 --normalizeUsing CPM --effectiveGenomeSize", eGenomeSize, ">", paste0(macs2dir, paste0("/", "group_", target_group[i]), "_bam2bw.log 2>&1")))
    }
}

cat(date(), paste("Pseudobulk analysis finished ...\n"))