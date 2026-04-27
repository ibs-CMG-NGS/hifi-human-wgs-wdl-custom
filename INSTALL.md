# 사전 요구사항 설치 가이드

Ubuntu 20.04 LTS 기준으로 작성되었습니다. 일부 단계는 root 권한(`sudo`)이 필요합니다.

## 목차

1. [git](#1-git)
2. [Conda (Miniconda)](#2-conda-miniconda)
3. [samtools](#3-samtools)
4. [Apptainer (Singularity)](#4-apptainer-singularity)
5. [NVIDIA 드라이버](#5-nvidia-드라이버)
6. [설치 확인](#6-설치-확인)

---

## 1. git

```bash
sudo apt update
sudo apt install -y git
git --version
```

---

## 2. Conda (Miniconda)

Conda는 Python 환경과 패키지를 관리합니다. Anaconda가 이미 설치되어 있다면 이 단계를 건너뜁니다.

```bash
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh
bash Miniconda3-latest-Linux-x86_64.sh -b -p $HOME/miniconda3
rm Miniconda3-latest-Linux-x86_64.sh
```

설치 후 셸에 conda 명령을 등록합니다:

```bash
$HOME/miniconda3/bin/conda init bash
source ~/.bashrc
conda --version
```

---

## 3. samtools

GRCm39 레퍼런스 FASTA 인덱스 생성(`samtools faidx`)에 필요합니다. conda 환경 생성 전에도 쓸 수 있도록 시스템에 설치합니다.

```bash
sudo apt install -y samtools
samtools --version
```

> apt 패키지가 구버전(1.x)일 경우 bioconda를 통해 최신 버전 설치를 권장합니다:
> ```bash
> conda install -c bioconda samtools
> ```

---

## 4. Apptainer (Singularity)

파이프라인의 모든 도구는 Apptainer 컨테이너 안에서 실행됩니다. Ubuntu 공식 PPA를 통해 설치합니다.

```bash
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt update
sudo apt install -y apptainer
apptainer --version
```

> PPA 사용이 불가능한 환경에서는 GitHub 릴리즈 페이지에서 `.deb` 파일을 직접 설치합니다:
> ```bash
> wget https://github.com/apptainer/apptainer/releases/download/v1.4.5/apptainer_1.4.5_amd64.deb
> sudo apt install ./apptainer_1.4.5_amd64.deb
> ```

Singularity 명령으로도 호출할 수 있도록 symlink가 필요합니다 (설치 시 자동 생성됨):

```bash
which singularity   # /usr/bin/singularity 로 나와야 함
```

---

## 5. NVIDIA 드라이버

GPU를 사용할 경우(DeepVariant 가속)에만 필요합니다. CPU만으로도 파이프라인 실행은 가능합니다.

### 현재 드라이버 상태 확인

```bash
nvidia-smi
```

출력이 나오면 드라이버가 이미 설치된 것입니다. 이 단계를 건너뜁니다.

### 드라이버 설치 (미설치 시)

Ubuntu가 GPU 모델에 맞는 드라이버를 자동으로 선택하게 합니다:

```bash
sudo apt update
sudo ubuntu-drivers autoinstall
sudo reboot
```

재부팅 후 확인:

```bash
nvidia-smi
```

### 특정 버전 지정 설치

자동 선택 대신 버전을 직접 지정하려면:

```bash
# 사용 가능한 드라이버 목록 확인
ubuntu-drivers devices

# 특정 버전 설치 (이 서버는 535 사용)
sudo apt install nvidia-driver-535
sudo reboot
```

### 이 서버의 GPU 구성

| 항목 | 내용 |
|------|------|
| GPU 모델 | NVIDIA GeForce RTX 2080 Ti × 2 |
| 드라이버 버전 | 535.183.01 |
| **주의** | GPU 0번 고장 — `CUDA_VISIBLE_DEVICES=1`로 GPU 1번만 사용 |

---

## 6. 설치 확인

모든 도구가 정상 설치되었는지 한 번에 확인합니다:

```bash
echo "=== git ===" && git --version
echo "=== conda ===" && conda --version
echo "=== samtools ===" && samtools --version | head -1
echo "=== apptainer ===" && apptainer --version
echo "=== singularity (symlink) ===" && singularity --version
echo "=== nvidia-smi ===" && nvidia-smi --query-gpu=name,driver_version --format=csv,noheader 2>/dev/null || echo "GPU 없음 또는 드라이버 미설치"
```

모두 확인되면 [SETUP.md](./SETUP.md)로 돌아가 파이프라인 환경을 구축합니다.
