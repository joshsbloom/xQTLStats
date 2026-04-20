# take the larger vcf,  genotype calls, and a subset of parents 
# extracts segregating sites 
# return an alphaSimR founder population
createFounderPop=function(vcf, gt, p.names, X.only=F, X.drop=T, gmap) { 
    gt.sub=gt[,colnames(gt) %in% p.names]

    #monomorphic=apply(gt.sub, 1, function(x) all.equal(x))
    #monomorphic sites 
    #faster to do this with math
    rSg=rowSums(gt.sub)
    #sites with hets 
    sum(is.na(rSg))
    #sites all ref
    sum(rSg==0, na.rm=T)
    #sites all alt
    sum(rSg==length(p.names), na.rm=T)
    #sites mito
    sum(grepl('MtDNA', rownames(gt.sub)))

    bad.sites= is.na(rSg) | rSg==0 | rSg==length(p.names)  | grepl('MtDNA', rownames(gt.sub))
    if(X.only) {  bad.sites = bad.sites | !(grepl('X_', rownames(gt.sub))) }
    if(X.drop) {  bad.sites = bad.sites | (grepl('X_', rownames(gt.sub))) }
    gt.sub=gt.sub[-which(bad.sites),]
    vcf.cross=vcf[match(rownames(gt.sub), rownames(gt)), samples=colnames(gt.sub)]
    #generate sample ID
    vcf.cross=vcfR::addID(vcf.cross)


    uchrU=unique(vcfR::getCHROM(vcf.cross))

    imputed.positions=jitterGmapVector(getGmapPositions(vcf.cross, gmap[uchrU], uchrU)) 
    

    #genetic map positions must be in Morgans
    genMap=data.frame(markerName=paste0(vcfR::getCHROM(vcf.cross),'_',vcfR::getPOS(vcf.cross)), 
                      chromosome=vcfR::getCHROM(vcf.cross), position=unlist(imputed.positions)/100)

    teg.GT=t(gt.sub)
    #recode
    teg.GT[teg.GT==0]=-1
    colnames(teg.GT)=paste0(vcfR::getCHROM(vcf.cross),'_',vcfR::getPOS(vcf.cross))
    ped=data.frame(id=rownames(teg.GT), mother=rep(0, nrow(teg.GT)), father=rep(0,nrow(teg.GT)) ) #c(0,0), father=c(0,0))
    return(AlphaSimR::importInbredGeno(geno=teg.GT, genMap=genMap, ped=ped))
}


#' Phase reference and alt counts given the ref/alt calls for one of the parents, using AlphaSimR objects
#'
#' @param df data.frame output from getBiallelicCounts() should contain the columns: ID,ref,alt
#' @param p1.name name of parent to be designated parent 1 in vcf.cross (name must exist in vcf.cross) 
#' @param founderPop #### 
#' @return data.frame of variant ID, ref count, and alt count phased  
#' @export
phaseBiparental=function(df, p1.name, founderPop, genMap){

        gID=genMap$id
        p1.ref=AlphaSimR::pullMarkerGeno(founderPop, genMap$id)[p1.name,]==0
        p1.ref=p1.ref[gID]
        vname=names(p1.ref)

        p1=c(df$ref[p1.ref], df$alt[!p1.ref])
        vscramb=c(vname[p1.ref], vname[!p1.ref])
        names(p1)=vscramb
        p1=p1[vname]

        p2=c(df$ref[!p1.ref], df$alt[p1.ref])
        vscramb=c(vname[!p1.ref], vname[p1.ref])
        names(p2)=vscramb
        p2=p2[vname]
        if(!is.null(df$expected)) {
            expected.phased=ifelse(p1.ref, df$expected, 1-df$expected)
            df$expected.phased=expected.phased
         }
         df$p1=p1
         df$p2=p2
    return(df)
}

#' Subset GATK table for known segregating variants and phase given parental genotypes
#'
#' @param sample_dir sample directory 
#' @param allele_counts name of file
#' @param p.names expected parents
#' @param vcf vcfR object
#' @param gt genotype calls from vcfR object
#' @param gmap genetic map
#' @return data.frame with original and phased counts
#' @export
makeCountTablesGATK=function(sample_dir, allele_counts, p.names, vcf, gt, gmap) {
	founderPop = createFounderPop(vcf,gt, p.names,X.only=F, X.drop=F, gmap)
	genMap=AlphaSimR::getGenMap(founderPop)
	scounts <- readr::read_tsv(stringr::str_c(sample_dir, "/", allele_counts))

	scounts.sub=scounts[paste0(scounts$contig, '_', scounts$position) %in% genMap$id,]
        scounts=data.frame(id=paste0(scounts.sub$contig, '_', scounts.sub$position),ref=scounts.sub$refCount, alt=scounts.sub$altCount)
        scounts=dplyr::left_join(genMap, scounts, by='id')

        scounts$ref[is.na(scounts$ref)]=0
        scounts$alt[is.na(scounts$alt)]=0
        names(scounts)[1]='ID'

	countdf=phaseBiparental(scounts, p.names[1], founderPop, genMap)

        #note, we need a better structure for keeping track of which parent is which
        attr(countdf, 'p1')=p.names[1]
        attr(countdf, 'p2')=p.names[2]

	return(countdf)

}


#' Subset GATK table for known segregating variants and phase given parental genotypes,
#' using a qs2 genotype object (Z_1011_filtered.qs) instead of vcfR/gt objects.
#'
#' @param sample_dir sample directory containing GATK ASEReadCounter output
#' @param allele_counts filename of the allele count table (e.g. "Sample001.txt")
#' @param p.names character vector of length 2: c(parent1, parent2). Parent1 is used as
#'   the phasing reference. Names must exist as rows in Z.
#' @param Z.obj list read from a qs2 genotype object (e.g. Z_1011_filtered.qs).
#'   Must contain: Z (dgCMatrix, strains x variants with 0=ref, 2=alt),
#'   variant_ids (character vector, format "chrI_27908:C/T"), strain_ids, chr.list.
#' @param gmap genetic map object (list of data.frames with ppos and map columns)
#' @param uchr vector of chromosome names in desired order (default: yeast chrI-XVI)
#' @return data.frame with original and phased counts
#' @export
makeCountTablesGATK2=function(sample_dir, allele_counts, p.names, Z.obj, gmap,
                              uchr=paste0('chr', as.character(as.roman(1:16)))) {

    # extract parent genotype vectors from the sparse Z matrix (0=ref, 2=alt, 1=het)
    p1.geno = Z.obj$Z[p.names[1], ]
    p2.geno = Z.obj$Z[p.names[2], ]

    # strip allele info from variant_ids: "chrI_27908:C/T" -> "chrI_27908"
    vid.short = sub(':.*$', '', Z.obj$variant_ids)

    # find segregating biallelic sites: both parents homozygous (0 or 2) and different
    seg = (p1.geno != p2.geno) & (p1.geno != 1) & (p2.geno != 1)
    seg.ids = vid.short[seg]
    p1.geno.seg = p1.geno[seg]

    # parse chromosome and position from the short variant IDs
    seg.chr = sub('_[0-9]+$', '', seg.ids)
    seg.pos = as.integer(sub('^.*_', '', seg.ids))

    # split positions by chromosome, impute genetic map positions
    p.by.chr = split(seg.pos, factor(seg.chr, levels = uchr))
    # only use chromosomes present in both gmap and data
    uchr.use = intersect(uchr, names(p.by.chr)[sapply(p.by.chr, length) > 0])
    p.by.chr = p.by.chr[uchr.use]

    imputed.positions = mapply(
        function(x, y) { approxfun(y$ppos, y$map, rule = 2)(x) },
        x = p.by.chr, y = gmap[uchr.use],
        SIMPLIFY = FALSE
    )
    imputed.positions = jitterGmapVector(imputed.positions)

    # build genMap-like data.frame (id, chromosome, position in Morgans)
    genMap = data.frame(
        id         = unlist(mapply(function(ids, chr) ids, ids = split(seg.ids, factor(seg.chr, levels = uchr.use)), chr = uchr.use, SIMPLIFY = FALSE)),
        chromosome = rep(uchr.use, sapply(p.by.chr, length)),
        position   = unlist(imputed.positions) / 100,
        stringsAsFactors = FALSE
    )

    # read GATK ASEReadCounter output and match to segregating sites
    scounts = readr::read_tsv(stringr::str_c(sample_dir, "/", allele_counts),
                              show_col_types = FALSE)
    scounts$id = paste0(scounts$contig, '_', scounts$position)
    scounts.sub = scounts[scounts$id %in% genMap$id, ]
    scounts = data.frame(id = scounts.sub$id,
                         ref = scounts.sub$refCount,
                         alt = scounts.sub$altCount)
    scounts = dplyr::left_join(genMap, scounts, by = 'id')
    scounts$ref[is.na(scounts$ref)] = 0
    scounts$alt[is.na(scounts$alt)] = 0
    names(scounts)[1] = 'ID'

    # phase: p1.ref is TRUE where parent1 carries the reference allele (geno == 0)
    p1.ref.lookup = (p1.geno.seg == 0)
    names(p1.ref.lookup) = seg.ids
    p1.ref = p1.ref.lookup[scounts$ID]
    vname = names(p1.ref)

    p1 = c(scounts$ref[p1.ref], scounts$alt[!p1.ref])
    vscramb = c(vname[p1.ref], vname[!p1.ref])
    names(p1) = vscramb
    p1 = p1[vname]

    p2 = c(scounts$ref[!p1.ref], scounts$alt[p1.ref])
    vscramb = c(vname[!p1.ref], vname[p1.ref])
    names(p2) = vscramb
    p2 = p2[vname]

    scounts$p1 = p1
    scounts$p2 = p2

    attr(scounts, 'p1') = p.names[1]
    attr(scounts, 'p2') = p.names[2]

    return(scounts)
}
