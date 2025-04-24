# 25년도 04월 24일 실습에 사용된 코드 정리본입니다.
## 결과 해석 말고, 코드 사용방법 등에 대해서는 도움을 드릴 수 있으니 11015khw@gmail.com으로 메일 넣어주시면 도와드릴 수 있는 부분은 도와드리겠습니다.

## Rtools설치가 필요할 수 있음.
## 현재 r version에 맞는 툴즈 설치하기
## https://cran.r-project.org/bin/windows/Rtools/

# install.packages("vcfR")
# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT", force=TRUE)
# install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# install.packages("ggpubr")

# ─────────────────────────────────────────────────────────────────────────────
# 1. 작업 환경 세팅 -------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# 현재 R 세션의 작업 디렉토리 확인
getwd()

# 프로젝트 루트 폴더 경로 지정
wk_dir <- "C:/Users/11015/Downloads/GAPIT_practice/GAPIT_practice/"

# 작업 디렉토리를 프로젝트 폴더로 변경
setwd(wk_dir)

# R 패키지 설치 (필요 시 한 번만 실행)
# devtools::install_github("jiabowang/GAPIT", force=TRUE)   # GAPIT 설치
# install.packages("vcfR")                                  # VCF 파일 처리용
# source("http://zzlab.net/GAPIT/gapit_functions.txt")      # GAPIT 함수 원격 로드
# install.packages("BiocManager")                           # Bioconductor 관리자
# BiocManager::install("SNPRelate")                         # SNPRelate 설치

# R 버전 정보 출력 (환경 확인용)
R.version


# ─────────────────────────────────────────────────────────────────────────────
# 2. VCF 파일 읽기 및 유전자형 매트릭스 추출 -----------------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(vcfR)

# VCF 포맷 파일(.vcf) 읽기
vcf <- read.vcfR(paste0(wk_dir, "data/mdp_genotype.vcf"))

# VCF 객체 자체 정보를 출력 (메타데이터, 샘플 수 등 확인)
vcf

# VCF 내용 미리보기: 상위 몇 줄만
head(vcf)

# GT(Genotype) 정보를 숫자형 매트릭스로 추출
# element="GT" : 유전자형, as.numeric=TRUE 로 0/1/2 형태로 변환
genotype_matrix <- extract.gt(vcf, element="GT", as.numeric=TRUE)

# 추출된 매트릭스 일부만 보기 (첫 5개 샘플 × 첫 5개 SNP)
head(genotype_matrix)[,1:5]

# 양대립형(biallelic) SNP만 필터링하여 확인
vcf[which(is.biallelic(vcf))]



# ─────────────────────────────────────────────────────────────────────────────
# 3. HapMap → 수치형 매트릭스 변환 (GAPIT) --------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(GAPIT)

# 표현형 데이터 읽기 (Taxa별 trait 값)
trait <- read.table(paste0(wk_dir, "data/mdp_traits.txt"), header=TRUE)

# HapMap 포맷의 유전자형 파일 읽기 (header=FALSE)
hapmap <- read.table(paste0(wk_dir, "data/mdp_genotype.hmp.txt"), header=FALSE)

# GAPIT::GAPIT() 호출 → 수치형 유전자형(GD), SNP 메타(GM) 생성
#   output.numerical=TRUE : HapMap → 0/1/2 matrix 변환
#   kinship.algorithm="Zhang" (기본) 방식으로 kinship 계산
genotype_matrix <- GAPIT(G = hapmap, output.numerical = TRUE)

# VanRaden 방식 kinship으로 바꿔서 다시 실행
genotype_matrix <- GAPIT(
  G = hapmap,
  output.numerical = TRUE,
  kinship.algorithm = "VanRaden"
)

# 생성된 결과 객체의 주요 요소 확인
names(genotype_matrix)
head(genotype_matrix$GM)        # SNP 메타데이터 (염색체·위치 등)
dim(genotype_matrix$GD[,-1])    # GD: 수치형 유전자형 (첫 열 taxa 제외)



# ─────────────────────────────────────────────────────────────────────────────
# 4. GDS 파일 생성 → LD pruning 및 PCA 준비 (SNPRelate) ------------------------
# ─────────────────────────────────────────────────────────────────────────────

library(SNPRelate)

# snpgdsCreateGeno(): HapMap → GDS 포맷으로 변환
snpgdsCreateGeno(
  gds.fn          = "output.gds",                                # 출력 파일명
  genmat          = as.matrix(genotype_matrix$GD[,-1]),         # 유전자형 matrix
  sample.id       = genotype_matrix$GD$taxa,                    # 샘플 ID
  snp.id          = genotype_matrix$GM$SNP,                     # SNP ID
  snp.rs.id       = genotype_matrix$GM$SNP,                     # rsID로도 동일 지정
  snp.chromosome  = genotype_matrix$GM$Chromosome,               # 염색체 정보
  snp.position    = genotype_matrix$GM$Position,                 # SNP 위치 (bp)
  snp.allele      = genotype_matrix$G$alleles,                  # 대립유전자 문자열
  snpfirstdim     = FALSE                                        # 행=샘플, 열=SNP
)

# 생성된 GDS 파일 열기
genofile <- snpgdsOpen("output.gds")

# LD pruning: 상관관계 기반으로 대표 SNP만 선별
snpset <- snpgdsLDpruning(
  gdsobj        = genofile,
  autosome.only = FALSE,    # 식물 유전체엔 성/비성 염색체 구분 없음
  remove.monosnp= TRUE,     # 변이 없는 SNP 제거
  maf           = 0.05,     # MAF 필터링
  missing.rate  = 0.1,      # 결측 허용률
  method        = "r",      # r² 기준
  slide.max.bp  = 10000L,   # 10kb 윈도우
  ld.threshold  = 0.2,      # r² 임계값
  start.pos     = "first",  # SNP 순서대로 진행
  num.thread    = 4L,       # 병렬 처리
  verbose       = TRUE
)

# 리스트 → 벡터 형태로 변환
snpset.id <- unlist(snpset)



# ─────────────────────────────────────────────────────────────────────────────
# 5. PCA 수행 ----------------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# snpgdsPCA(): 선별된 SNP로 PCA
pca <- snpgdsPCA(
  gdsobj        = genofile,
  snp.id        = snpset.id,
  autosome.only = FALSE,
  remove.monosnp= TRUE,
  maf           = 0.05,
  missing.rate  = 0.1,
  algorithm     = "exact",   # 정확 계산
  eigen.cnt     = 10L,       # 10개 PC 계산
  num.thread    = 4L,
  bayesian      = FALSE,
  verbose       = TRUE
)

# 각 PC가 설명하는 분산 비율(%) 계산
pc_percent <- pca$varprop[1:10] * 100

# Scree plot: PC별 분산비 시각화
data.frame(
  PC       = factor(paste0("PC", 1:10), levels = paste0("PC", 1:10)),
  variance = pc_percent
) |>
  ggplot2::ggplot(ggplot2::aes(PC, variance)) +
  ggplot2::geom_bar(stat="identity", fill="skyblue") +
  ggplot2::geom_line(ggplot2::aes(group=1), color="blue", linewidth=1) +
  ggplot2::geom_point(color="darkblue", size=2) +
  ggplot2::labs(
    x     = "Principal Component",
    y     = "Explained Variance (%)",
    title = "Scree Plot"
  ) +
  ggplot2::theme_minimal()



# ─────────────────────────────────────────────────────────────────────────────
# 6. 군집 분석: NJ 트리 생성 --------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# GAPIT 내장 기능으로 NJ tree 기반 그룹화 (3개 군집)
group_analysis <- GAPIT(
  G               = hapmap,
  output.numerical= TRUE,
  NJtree.group    = 3
)

# 생성된 군집 결과 파일 읽기
group <- read.table(
  paste0(wk_dir, "interim_results/GAPIT.Genotype.Kin_NJtree_compress_z.txt"),
  header = TRUE
)



# ─────────────────────────────────────────────────────────────────────────────
# 7. PCA 결과 + 군집 정보 통합 & 시각화 ----------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

# PCA 결과와 NJ 군집 번호를 하나의 데이터프레임으로 통합
pca_table <- data.frame(
  sample = pca$sample.id,
  group  = factor(t(group) %*% matrix(c(1,2,3), ncol=1)),
  PC1     = pca$eigenvect[,1],  # PC1 로딩
  PC2     = pca$eigenvect[,2],  # PC2 로딩
  PC3     = pca$eigenvect[,3]   # PC3 로딩
)

# 2D 산점도: PC1 vs PC2
p1 <- ggplot2::ggplot(pca_table, ggplot2::aes(PC1, PC2, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x=paste0("PC1 (", round(pc_percent[1],2), "%)"),
    y=paste0("PC2 (", round(pc_percent[2],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position="none")

# 2D 산점도: PC2 vs PC3
p2 <- ggplot2::ggplot(pca_table, ggplot2::aes(PC2, PC3, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x=paste0("PC2 (", round(pc_percent[2],2), "%)"),
    y=paste0("PC3 (", round(pc_percent[3],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position="none")

# 2D 산점도: PC1 vs PC3
p3 <- ggplot2::ggplot(pca_table, ggplot2::aes(PC1, PC3, color=group)) +
  ggplot2::geom_point(size=3, alpha=0.8) +
  ggplot2::labs(
    x=paste0("PC1 (", round(pc_percent[1],2), "%)"),
    y=paste0("PC3 (", round(pc_percent[3],2), "%)")
  ) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.position="none")

# 세 개 플롯 한 화면에 배치
ggpubr::ggarrange(p1, p2, p3, ncol=3)

# 3D 산점도 (Plotly)
plotly::plot_ly(
  pca_table,
  x = ~PC1, y = ~PC2, z = ~PC3,
  color = ~group,
  type = "scatter3d",
  mode = "markers",
  marker = list(size = 5)
)



# ─────────────────────────────────────────────────────────────────────────────
# 8. GWAS 실행 (GLM, MLM) ------------------------------------------------------
# ─────────────────────────────────────────────────────────────────────────────

gwas_result <- GAPIT(
  Y     = trait[, c(1,2)],      # 표현형: Taxa 와 trait 값
  G     = hapmap,               # HapMap 유전자형
  CV    = data.frame(           # 공변량: 상위 3개 PC
    Taxa = pca$sample.id,
    PC1  = pca$eigenvect[,1],
    PC2  = pca$eigenvect[,2],
    PC3  = pca$eigenvect[,3]
  ),
  model = c("GLM", "MLM")        # 사용할 모델: GLM, MLM
)

