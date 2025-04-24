# RDA_preparation
Repository for 2025_04_24 RDA lecture.
* 해당 repository는 결과 해석에 대한 내용은 포함하지 않습니다.
* 결과 해석 부분은 아래 related links의 3번에서 확인 부탁드립니다.

아래 qr 링크: [https://url.kr/evj1tb]

![QrCode-24-04-2025_11-00-22](https://github.com/user-attachments/assets/30c57d4b-bb48-4e29-9a5f-f2f8da5c985b)


## related links   
### Original Author's github   
All the following lecture materials were written by Professor Seungyong Hwang.    
* [https://github.com/vic-dragon]    
* [https://vic-dragon.github.io/]

### 1. lecture PDF  
[https://url.kr/5d7wp8]

### 2. Dataset for this lecture
[https://url.kr/wrs6lb]

### 3. Codes
[https://url.kr/tz8gs7]

```
NOTE: The three links(1. ~ 3.) above will expire on April 21, 2028.
```

## Install Libraries
``` r
## Rtools설치가 필요할 수 있음.
## 현재 r version에 맞는 툴즈 설치하기
## https://cran.r-project.org/bin/windows/Rtools/

# install.packages("vcfR")
# install.packages("devtools")
# devtools::install_github("jiabowang/GAPIT", force=TRUE)
# install.packages("BiocManager")
# BiocManager::install("SNPRelate")
# install.packages("ggpubr")
```
NOTE: You have to Use R version =< 4.4.1 to use GAPIT library   

집(개인 수준으로)에서 하실 때에는, 
* version을 확인 후 4.4.1 이하의 R 버전을 깔았는지 확인을 한번 해야 합니다.
저 같은 경우는 R 4.3.3 또는 4.3.1을 사용하고 있습니다.

아래 링크에서 우선 해당 버전을 설치해주시고
* [https://cran.r-project.org/bin/windows/base/old/]

아래 링크에서 RTools 4.3를 설치해주시면 됩니다.
* [https://cran.r-project.org/bin/windows/Rtools/]


