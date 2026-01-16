# 설치 가이드

## 1. Conda 환경 설정

### 환경 생성
```bash
cd /home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL

# conda 환경 생성
conda env create -f environment.yml

# 환경 활성화
conda activate hifi-human-wgs

# 설치 확인
miniwdl --version
```

또는 자동 설정 스크립트 사용:
```bash
chmod +x setup_environment.sh
./setup_environment.sh
```

## 2. Singularity 설치 (필수)

Singularity는 conda로 설치할 수 없으므로 시스템에 별도로 설치해야 합니다.

### Ubuntu/Debian 시스템에서 설치

#### 방법 1: 패키지 관리자 사용 (권장)
```bash
# SingularityCE 최신 버전 설치
sudo apt-get update
sudo apt-get install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt-get update
sudo apt-get install -y apptainer
```

Apptainer는 Singularity의 후속 프로젝트이며 완벽히 호환됩니다.

#### 방법 2: 소스에서 빌드
```bash
# 의존성 패키지 설치
sudo apt-get update
sudo apt-get install -y \
    build-essential \
    libseccomp-dev \
    pkg-config \
    squashfs-tools \
    cryptsetup \
    curl wget git

# Go 설치 (필요시)
export VERSION=1.21.0 OS=linux ARCH=amd64
wget https://dl.google.com/go/go$VERSION.$OS-$ARCH.tar.gz
sudo tar -C /usr/local -xzvf go$VERSION.$OS-$ARCH.tar.gz
export PATH=/usr/local/go/bin:$PATH

# SingularityCE 설치
export VERSION=4.0.2
wget https://github.com/sylabs/singularity/releases/download/v${VERSION}/singularity-ce-${VERSION}.tar.gz
tar -xzf singularity-ce-${VERSION}.tar.gz
cd singularity-ce-${VERSION}
./mconfig
make -C builddir
sudo make -C builddir install
```

### 설치 확인
```bash
singularity --version
# 또는
apptainer --version
```

## 3. miniwdl 설정

### 설정 파일 복사 및 수정
```bash
# 설정 디렉토리 생성
mkdir -p ./config

# 예제 설정 파일 복사
cp backends/hpc/miniwdl.cfg ./config/miniwdl.cfg

# 설정 파일 편집
nano ~/.config/miniwdl.cfg
```

### 주요 설정 항목

`~/.config/miniwdl.cfg` 파일에서 다음 항목들을 확인/수정하세요:

```ini
[file_io]
# map 파일 사용을 위해 필수
allow_any_input = true

[singularity]
# Singularity/Apptainer 실행 파일 경로
exe = ["singularity"]  # 또는 ["apptainer"]

[scheduler]
# SLURM 사용시
container_backend = slurm_singularity
```

## 4. 참조 데이터 다운로드

워크플로우 실행을 위해 참조 데이터가 필요합니다:

```bash
# 작업 디렉토리 생성
mkdir -p ~/hifi-wdl-resources
cd ~/hifi-wdl-resources

# 참조 데이터 다운로드 (약 50GB)
wget https://zenodo.org/record/17086906/files/hifi-wdl-resources-v3.1.0.tar

# 압축 해제
tar -xvf hifi-wdl-resources-v3.1.0.tar

# 경로 기록
echo "참조 데이터 경로: $(pwd)"
```

## 5. 입력 파일 설정

```bash
cd /home/ygkim/ngs_pipeline/HiFi-human-WGS-WDL

# HPC 입력 템플릿 복사
cp backends/hpc/singleton.hpc.inputs.json my_inputs.json

# 입력 파일 편집
nano my_inputs.json

# <local_path_prefix>를 실제 참조 데이터 경로로 변경
# 예: /home/ygkim/hifi-wdl-resources
```

## 6. 워크플로우 실행 테스트

```bash
# conda 환경 활성화
conda activate hifi-human-wgs

# 워크플로우 실행 (드라이런)
miniwdl run workflows/singleton.wdl --input my_inputs.json --verbose

# 실제 실행
miniwdl run workflows/singleton.wdl --input my_inputs.json
```

## 7. 문제 해결

### Singularity/Apptainer를 찾을 수 없는 경우
```bash
# 설치 확인
which singularity
which apptainer

# PATH에 추가 (필요시)
export PATH=/usr/local/bin:$PATH

# ~/.bashrc에 영구 추가
echo 'export PATH=/usr/local/bin:$PATH' >> ~/.bashrc
source ~/.bashrc
```

### miniwdl이 Singularity를 인식하지 못하는 경우
```bash
# ~/.config/miniwdl.cfg 확인
cat ~/.config/miniwdl.cfg

# singularity 섹션 확인
[singularity]
exe = ["singularity"]  # 또는 apptainer가 설치된 경우 ["apptainer"]
```

### 인터넷 접속이 제한된 컴퓨트 노드 사용시
```bash
# 로그인 노드에서 미리 이미지 캐시 준비
./scripts/populate_miniwdl_singularity_cache.sh
```

## 8. 추가 리소스

- [HPC 백엔드 문서](docs/backend-hpc.md)
- [miniwdl 문서](https://miniwdl.readthedocs.io/)
- [Singularity/Apptainer 문서](https://apptainer.org/docs/)
- [메인 README](README.md)

## 요약

1. ✅ Conda 환경 생성: `conda env create -f environment.yml`
2. ✅ Singularity 설치: `sudo apt-get install apptainer`
3. ✅ miniwdl 설정: `cp backends/hpc/miniwdl.cfg ~/.config/`
4. ✅ 참조 데이터 다운로드 및 압축 해제
5. ✅ 입력 파일 설정 및 경로 수정
6. ✅ 워크플로우 실행 테스트
