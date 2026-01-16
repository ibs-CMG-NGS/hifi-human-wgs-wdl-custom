# Git Repository Setup - Quick Start

## 초기 설정

### 1. 현재 상태 확인
```bash
git status
```

### 2. 개인 설정 파일 생성
실제 데이터 경로가 포함된 설정 파일은 git에 추적되지 않습니다:

```bash
# 예제 파일을 복사하여 개인 설정 생성
cp sample.inputs.json.example my_sample.inputs.json

# 개인 설정 파일 편집
nano my_sample.inputs.json
```

### 3. Git 저장소 초기화 (아직 안했다면)
```bash
git init
git add .
git commit -m "Initial commit: HiFi WGS pipeline setup"
```

### 4. 리모트 저장소 연결
```bash
git remote add origin <your-repo-url>
git push -u origin main
```

## 추적되지 않는 파일들 (`.gitignore`에 포함)

### 데이터 파일
- `data/` - 모든 시퀀싱 데이터
- `*.bam`, `*.bai` - 대용량 BAM 파일들
- `hifi-wdl-resources/` - 다운로드된 참조 파일들

### 실행 결과
- `20*/` - 타임스탬프가 있는 워크플로우 실행 디렉토리
- `call-*/` - 개별 작업 실행 결과
- `error.json`, `rerun` - 에러 및 재실행 파일

### 캐시
- `miniwdl_call_cache/`
- `miniwdl_singularity_cache/`
- `__pycache__/`

### 개인 설정
- `*.inputs.json` - 개인 입력 설정 파일
- `*.local.*` - 로컬 설정 오버라이드
- `.env` - 환경 변수

## 추적되는 중요 파일들

### 워크플로우 정의
- `workflows/*.wdl` - 워크플로우 정의
- `workflows/*.inputs.json` - 기본 템플릿

### 설정 템플릿
- `*.inputs.json.example` - 입력 파일 예제
- `GRCh38.*.template.tsv` - 참조 파일 템플릿
- `config/miniwdl.cfg` - 워크플로우 엔진 설정

### 환경 설정
- `environment.yml` - Conda 환경
- `requirements.txt` - Python 의존성
- `setup_environment.sh` - 환경 설정 스크립트

### 문서
- `README*.md` - 프로젝트 문서
- `docs/` - 상세 문서
- `CONFIG_MANAGEMENT.md` - 설정 관리 가이드

## 팀 협업 시 권장사항

1. **개인 설정 공유하지 않기**: `*.inputs.json` 파일은 개인 경로를 포함하므로 커밋하지 마세요.

2. **예제 파일 업데이트**: 새로운 설정 옵션을 추가할 때는 `*.example` 파일도 함께 업데이트하세요.

3. **문서화**: 로컬 환경 설정이 필요한 경우 README에 명확히 기록하세요.

4. **상대 경로 사용**: 가능한 경우 절대 경로 대신 상대 경로를 사용하세요.

## 자주 사용하는 명령어

```bash
# 변경사항 확인
git status

# 새 파일 추가 및 커밋
git add <file>
git commit -m "설명"

# 원격 저장소에 푸시
git push

# 최신 변경사항 가져오기
git pull

# 무시된 파일 확인 (디버깅용)
git status --ignored
```

## 문제 해결

### 실수로 대용량 파일을 커밋한 경우
```bash
# 마지막 커밋 취소 (아직 push 안한 경우)
git reset HEAD~1

# 이미 push한 경우 (주의: 히스토리 재작성)
git filter-branch --force --index-filter \
  "git rm --cached --ignore-unmatch <파일경로>" \
  --prune-empty --tag-name-filter cat -- --all
```

### .gitignore 규칙이 작동하지 않는 경우
```bash
# 이미 추적 중인 파일은 캐시에서 제거 필요
git rm -r --cached .
git add .
git commit -m "Update .gitignore"
```

## 추가 정보

자세한 설정 관리 방법은 `CONFIG_MANAGEMENT.md`를 참조하세요.
