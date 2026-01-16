# Conda 환경 설정 가이드

이 프로젝트는 HiFi-human-WGS-WDL 워크플로우를 실행하기 위한 conda 환경을 제공합니다.

## 환경 설정

### 방법 1: 자동 설정 스크립트 사용

```bash
# 스크립트에 실행 권한 부여
chmod +x setup_environment.sh

# 환경 설정 실행
./setup_environment.sh
```

### 방법 2: 수동 설정

```bash
# environment.yml 파일을 사용하여 conda 환경 생성
conda env create -f environment.yml

# 환경 활성화
conda activate hifi-human-wgs

# 설치 확인
miniwdl --version
```

### 방법 3: pip를 사용한 최소 설정

conda를 사용하지 않는 경우:

```bash
# Python 가상환경 생성 (선택사항)
python -m venv venv
source venv/bin/activate  # Windows의 경우: venv\Scripts\activate

# requirements.txt로 패키지 설치
pip install -r requirements.txt
```

## 포함된 도구

환경에는 다음 도구들이 포함되어 있습니다:

- **miniwdl** (>=1.9.0): WDL 워크플로우 실행 엔진
- **miniwdl-slurm**: HPC 백엔드를 위한 SLURM 플러그인
- **singularity-ce**: 컨테이너 런타임
- **cromwell-tools**: Cromwell 워크플로우 관리 도구

## 환경 사용

환경을 활성화한 후 워크플로우를 실행할 수 있습니다:

```bash
# 환경 활성화
conda activate hifi-human-wgs

# miniwdl을 사용한 워크플로우 실행
miniwdl run workflows/singleton.wdl --input backends/hpc/singleton.hpc.inputs.json

# 또는 family 워크플로우 실행
miniwdl run workflows/family.wdl --input backends/hpc/family.hpc.inputs.json
```

## 환경 관리

```bash
# 환경 업데이트
conda env update -n hifi-human-wgs -f environment.yml --prune

# 환경 제거
conda env remove -n hifi-human-wgs

# 설치된 패키지 목록 확인
conda list -n hifi-human-wgs
```

## 추가 설정

### miniwdl 설정

HPC 백엔드를 사용하는 경우, miniwdl 설정 파일을 복사하고 편집하세요:

```bash
# 설정 디렉토리 생성
mkdir -p ~/.config

# 예제 설정 파일 복사
cp backends/hpc/miniwdl.cfg ~/.config/miniwdl.cfg

# 설정 파일 편집 (SLURM 설정에 맞게 수정)
nano ~/.config/miniwdl.cfg
```

### Singularity 캐시 설정

인터넷 접속이 제한된 컴퓨트 노드를 사용하는 경우:

```bash
# Singularity 캐시 사전 준비
./scripts/populate_miniwdl_singularity_cache.sh
```

## 참조 데이터 다운로드

워크플로우 실행을 위해 참조 데이터를 다운로드하세요:

```bash
# 참조 데이터 다운로드
wget https://zenodo.org/record/17086906/files/hifi-wdl-resources-v3.1.0.tar

# 압축 해제
tar -xvf hifi-wdl-resources-v3.1.0.tar
```

## 문제 해결

### miniwdl 설치 오류

conda-forge 채널에서 miniwdl을 찾을 수 없는 경우:

```bash
conda activate hifi-human-wgs
pip install miniwdl>=1.9.0 miniwdl-slurm
```

### Singularity 설치 문제

시스템에 Singularity가 이미 설치되어 있는 경우, environment.yml에서 singularity-ce를 제거할 수 있습니다.

## 추가 정보

- [HPC 백엔드 문서](docs/backend-hpc.md)
- [Singleton 워크플로우 문서](docs/singleton.md)
- [Family 워크플로우 문서](docs/family.md)
- [메인 README](README.md)
