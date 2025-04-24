# 작업 환경 세팅 -------------------------------------------------------------
# 현재 작업 디렉토리 확인
getwd()

# 작업 디렉토리 경로 지정 (GAPIT 실습 폴더)
wk_dir <- "C:/Users/11015/Downloads/GAPIT_practice/GAPIT_practice/"
# 작업 디렉토리 변경
setwd(wk_dir)

# 필요한 패키지 설치 (한 번만 실행)
# devtools::install_github("jiabowang/GAPIT", force=TRUE)   # GAPIT 패키지 설치
# install.packages("vcfR")                                  # VCF 파일 파싱용 vcfR
# source("http://zzlab.net/GAPIT/gapit_functions.txt")      # GAPIT 함수 원격 소스 로드
# install.packages("BiocManager")                           # Bioconductor 관리자
# BiocManager::install("SNPRelate")                         # SNPRelate 패키지 설치


# VCF 파일 불러오기 -----------------------------------------------------------
# vcfR::read.vcfR() 로 VCF 포맷의 유전체 변이 정보 읽기
# file 경로: data/mdp_genotype.vcf
vcf <- vcfR::read.vcfR(paste0(wk_dir, "data/mdp_genotype.vcf"))

# VCF 객체 출력하여 정보 확인
vcf
# 상위 몇 줄 출력하여 데이터 구조 살펴보기
head(vcf)  # VCF 파일 프리뷰


# Genotype matrix 추출 --------------------------------------------------------
# extract.gt() 로 genotype(유전자형) 정보만 수치형 matrix 로 변환
genotype_matrix <- vcfR::extract.gt(vcf, element="GT", as.numeric=TRUE)
# element 옵션:
#  - "GT": genotype (유전자형)
#  - "DP": read depth (깊이)
#  - "GQ": genotype quality (품질)
# 상위 5개 종(sample) × 5개 SNP 확인
head(genotype_matrix)[, 1:5]

# Biallelic SNP만 필터링 (양대립형만)
vcf[which(vcfR::is.biallelic(vcf))]


# LD pruning 및 PCA 준비 ------------------------------------------------------
library(GAPIT)
# 표현형(trait) 데이터 읽기 (mdp_traits.txt: 종별 형질값)
trait <- read.table(paste0(wk_dir, "data/mdp_traits.txt"), header=TRUE)
# HapMap 형식의 유전자형 파일 읽기
hapmap <- read.table(paste0(wk_dir, "data/mdp_genotype.hmp.txt"), header=FALSE)

# GAPIT 으로 수치형(genotype matrix) 생성 (기본 kinship: Zhang)
genotype_matrix <- GAPIT::GAPIT(G = hapmap, output.numerical = TRUE)

# kinship 알고리즘을 VanRaden 방식으로 변경하여 재실행
genotype_matrix <- GAPIT::GAPIT(
  G = hapmap,
  output.numerical = TRUE,
  kinship.algorithm = 'VanRaden'
)

# 결과 객체 구조 확인
names(genotype_matrix)
# GM: SNP 메타데이터 (chromosome, position 등)
head(genotype_matrix$GM)
# GD: 수치형 유전자형 데이터 (첫 열은 taxa)
dim(genotype_matrix$GD[, -1])


# GDS 파일 생성 (SNPRelate) -------------------------------------------------
# snpgdsCreateGeno() 로 GDS 포맷 파일 생성
SNPRelate::snpgdsCreateGeno(
  gds.fn = 'output.gds',                                 # 출력 파일명
  genmat = as.matrix(genotype_matrix$GD[, -1]),          # genotype matrix
  sample.id = genotype_matrix$GD$taxa,                   # 샘플 ID
  snp.id = genotype_matrix$GM$SNP,                       # SNP ID
  snp.rs.id = genotype_matrix$GM$SNP,                    # rsID
  snp.chromosome = genotype_matrix$GM$Chromosome,        # 염색체 번호
  snp.position = genotype_matrix$GM$Position,            # SNP 위치(bp)
  snp.allele = genotype_matrix$G$alleles,                # 대립유전자 정보
  snpfirstdim = FALSE                                    # 행: sample, 열: SNP
)

# 생성한 GDS 파일 열기
genofile <- SNPRelate::snpgdsOpen(paste0(wk_dir, "data/output.gds"))


# LD pruning (상관관계 기반 SNP 축소) -----------------------------------------
snpset <- SNPRelate::snpgdsLDpruning(
  gdsobj         = genofile,
  autosome.only  = FALSE,    # 식물에는 성/비성 염색체 구분 없음
  remove.monosnp = TRUE,     # 변이가 없는 SNP 제거
  maf            = 0.05,     # minor allele frequency 기준
  missing.rate   = 0.1,      # 결측치 허용률
  method         = "r",      # 상관계수 r² 기준
  slide.max.bp   = 10000L,   # 10kb 윈도우
  ld.threshold   = 0.2,      # r² 임계값
  start.pos      = "first",   # 순차적 선택
  num.thread     = 4L,       # 병렬 처리
  verbose        = TRUE
)
# pruning 후 선택된 SNP ID 리스트
snpset.id <- unlist(snpset)


# PCA 분석 (주성분 분석) ------------------------------------------------------
pca <- SNPRelate::snpgdsPCA(
  gdsobj         = genofile,
  snp.id         = snpset.id,
  autosome.only  = FALSE,
  remove.monosnp = TRUE,
  maf            = 0.05,
  missing.rate   = 0.1,
  algorithm      = "exact",  # 정확 계산
  eigen.cnt      = 10L,      # 추출할 주성분 수
  num.thread     = 4L,
  bayesian       = FALSE,
  verbose        = TRUE
)

# 각 PC의 분산 비율(%) 계산
pc_percent <- pca$varprop[1:10] * 100
# Scree plot 그리기
data.frame(
  PC       = factor(paste0('PC', 1:length(pc_percent)), levels = paste0('PC', 1:length(pc_percent))),
  variance = pc_percent
) %>%
  ggplot2::ggplot(ggplot2::aes(x=PC, y=variance)) +
  ggplot2::geom_bar(stat='identity', fill='skyblue') +
  ggplot2::geom_line(ggplot2::aes(group=1), linewidth=1) +
  ggplot2::geom_point(size=2) +
  ggplot2::labs(
    x = 'Principal Component',
    y = 'Explained Variance (%)',
    title = 'Scree Plot'
  ) +
  ggplot2::theme_minimal()


# 군집 분석 (NJ tree) --------------------------------------------------------
group_analysis <- GAPIT::GAPIT(
  G             = hapmap,
  output.numerical = TRUE,
  NJtree.group  = 3       # 3개의 군집으로 NJ tree 생성
)

# NJ tree 결과 파일 읽기 (중간 결과)
group <- read.table(
  paste0(wk_dir, "interim_results/GAPIT.Genotype.Kin_NJtree_compress_z.txt"),
  header = TRUE
)


# PCA 결과와 그룹 정보 통합 -----------------------------------------------
pca_table <- data.frame(
  sample = pca$sample.id,
  group  = factor(t(group) %*% matrix(c(1,2,3), ncol=1)),  # 그룹 번호 할당
  PC1    = pca$eigenvect[,1],  # PC1 값
  PC2    = pca$eigenvect[,2],  # PC2 값
  PC3    = pca$eigenvect[,3]   # PC3 값
)

# 2D scatter plots (PC1 vs PC2, PC2 vs PC3, PC1 vs PC3)
p1 <- ggplot2::ggplot(pca_table, ggplot2::aes(x=PC1, y=PC2, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x = paste0("PC1 (", round(pc_percent[1],2), "%)"),
    y = paste0("PC2 (", round(pc_percent[2],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position='none')

p2 <- ggplot2::ggplot(pca_table, ggplot2::aes(x=PC2, y=PC3, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x = paste0("PC2 (", round(pc_percent[2],2), "%)"),
    y = paste0("PC3 (", round(pc_percent[3],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position='none')

p3 <- ggplot2::ggplot(pca_table, ggplot2::aes(x=PC1, y=PC3, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x = paste0("PC1 (", round(pc_percent[1],2), "%)"),
    y = paste0("PC3 (", round(pc_percent[3],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position='none')

# 세 개의 플롯을 한 화면에 배치
ggpubr::ggarrange(p1, p2, p3, ncol=3)


# 3D PCA plot (Plotly) -------------------------------------------------------
plotly::plot_ly(
  pca_table,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~group,
  type = "scatter3d", mode = "markers",
  marker = list(size = 5)
)


# GWAS 수행 (Genome-Wide Association Study) --------------------------------
gwas_result <- GAPIT::GAPIT(
  Y     = trait[, c(1,2)],    # 표현형 데이터: 첫 두 열 (Taxa, trait 값)
  G     = hapmap,             # 유전자형 HapMap
  CV    = data.frame(         # 공변량: PCA 상위 3개 PC
    Taxa = pca$sample.id,
    PC1  = pca$eigenvect[,1],
    PC2  = pca$eigenvect[,2],
    PC3  = pca$eigenvect[,3]
  ),
  model = c('GLM','MLM')       # 사용할 모델: GLM, MLM 등
)
