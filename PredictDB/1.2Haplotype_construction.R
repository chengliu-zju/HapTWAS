# 定义单倍型
ld <- 0.8
plink<-'plink'
plink2<-'plink2'
library('GHap')
library("data.table")
for (i in 1:22) {
  #用plink计算LD定义haplotype block
  system(paste0(plink,' --vcf ./data/SNP/chr',i,'.vcf.gz --blocks no-pheno-req --blocks-max-kb 200 --blocks-strong-lowci ',ld,' --double-id --out chr',i))
  #读入上一步生成的.blocks.det文件并处理为GHap需要的格式（"BLOCK","CHR","BP1","BP2","SIZE","NSNP"）
  blocks.ld <- read.table(paste0("chr",i,".blocks.det"), header=T)[,1:5]
  blocks.ld <- cbind(blocks.ld[,1],blocks.ld)
  blocks.ld[,5] <- blocks.ld[,5]*1000
  rownames(blocks.ld) <- 1:nrow(blocks.ld)
  blocks.ld[,1] <- paste0('CHR',blocks.ld[,1],'_B',rownames(blocks.ld))
  colnames(blocks.ld) <-c ("BLOCK","CHR","BP1","BP2","SIZE","NSNP")
  #plink2生成sample and haps 两个文件用于下一步输入
  system(paste0(plink2,' --vcf ./data/SNP/chr',i,'.vcf.gz --export haps --out chr',i))
  #ghap.oxford2phase函数生成GHap三个输入文件(samples包括population and ID即fam文件，markers包括Chromosome,Marker ID,Position(bp),Ref Allele(A0) and Alt Allele(A1)五列即bim文件，phase有m行（标记）2n列（样本)即ped文件基因型用0/1编码）
  ghap.oxford2phase(input.files = paste0("chr",i), out.file = paste0("chr",i),ncores=10)
  #压缩文件，输入上一步的三个文件构建GHap phase对象（.phaseb文件）
  ghap.compress(input.file = paste0("chr",i), out.file = paste0("chr",i),ncores = 10)
  #读入.phaseb文件
  phase <- ghap.loadphase(paste0("chr",i))
  #构建基于haplotype的基因型矩阵,输出hapgenotypesb,hapalleles and hapsamples文件（binary二进制输出，ncores线程）
  ghap.haplotyping(phase, blocks.ld, outfile = paste0("chr",i), binary = T, ncores = 10)
  #读入haplotype基因型
  haplo <- ghap.loadhaplo(paste0("chr",i))
  #数据量大时很耗时，将haplotype输出为plink格式(BED BIM FAM)
  ghap.hap2plink(haplo, outfile = paste0("chr",i))
  #将bed bim fam 转码为样本对单体型012编码的raw文件（类似于一对姐妹染色单体的aa/Aa/AA）
  system(paste0(plink," --bfile chr",i," --recodeA --out chr",i))
}

sw <- 20
plink<-'plink'
plink2<-'plink2'
library('GHap')
library("data.table")
for (i in 1:22) {
  #plink2生成sample and haps 两个文件用于下一步输入
  system(paste0(plink2,' --vcf ./data/SNP/chr',i,'.vcf.gz --export haps --out chr',i))
  #ghap.oxford2phase函数生成GHap三个输入文件(samples包括population and ID即fam文件，markers包括Chromosome,Marker ID,Position(bp),Ref Allele(A0) and Alt Allele(A1)五列即bim文件，phase有m行（标记）2n列（样本)即ped文件基因型用0/1编码）
  ghap.oxford2phase(input.files = paste0("chr",i), out.file = paste0("chr",i), ncores=10)
  #压缩文件，输入上一步的三个文件构建GHap phase对象（.phaseb文件）
  ghap.compress(input.file = paste0("chr",i), out.file = paste0("chr",i), ncores = 10)
  #读入.phaseb文件
  phase <- ghap.loadphase(paste0("chr",i))
  # 用marker个数进行sliding
  blocks.mkr <- ghap.blockgen(phase, windowsize = sw, slide = sw, unit = "marker")
  # 构建基于haplotype的基因型矩阵,输出hapgenotypesb,hapalleles and hapsamples文件（binary二进制输出，ncores线程）
  ghap.haplotyping(phase, blocks.mkr, outfile = paste0("chr",i), binary = T, ncores = 10)
  #读入haplotype基因型
  haplo <- ghap.loadhaplo(paste0("chr",i))
  #数据量大时很耗时，将haplotype输出为plink格式(BED BIM FAM)
  ghap.hap2plink(haplo, outfile = paste0("chr",i))
  #将bed bim fam 转码为样本对单体型012编码的raw文件（类似于一对姐妹染色单体的aa/Aa/AA）
  system(paste0(plink," --bfile chr",i," --recodeA --out chr",i))
}

