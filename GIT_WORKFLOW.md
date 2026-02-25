# Git 워크플로우

## .gitignore 원칙

### Git에 포함되지 않는 것 (개인/대용량)

- `*.inputs.json` — 개인 데이터 경로 포함
- `*.local.*` — 로컬 설정 오버라이드
- `hifi-wdl-resources/` — 대용량 레퍼런스 데이터
- `data/`, `*.bam`, `*.bai` — 시퀀싱 데이터
- `miniwdl_call_cache/`, `miniwdl_singularity_cache/` — 실행 캐시
- `batch_results/`, `_2*/` — 분석 결과

### Git에 포함되는 것 (공유/템플릿)

- `workflows/*.wdl` — 워크플로우 정의
- `config/miniwdl.cfg` — 기본 설정 (local.cfg는 제외)
- `*.inputs.json.example` — 입력 파일 템플릿
- `*.ref_map.*.template.tsv` — 레퍼런스 맵 템플릿
- `environment.yml`, `requirements.txt` — 환경 명세
- `*.md` — 문서
- `*.sh` — 스크립트

## 로컬 변경사항 커밋 (서버에서)

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom

# 변경 파일 확인
git status

# 특정 파일 추가
git add config/miniwdl.cfg
git add SETUP.md BATCH_GUIDE.md SERVER_GUIDE.md GIT_WORKFLOW.md TROUBLESHOOTING.md

# 커밋
git commit -m "설명 메시지"

# 원격 저장소에 푸시
git push origin main
```

## Windows/VS Code → 서버 동기화

### Windows에서 push:

```powershell
cd \\wsl.localhost\Ubuntu\home\ygkim\ngs-pipeline\HiFi-human-WGS-WDL
git add *.sh *.md config/
git commit -m "Update scripts and docs"
git push origin main
```

### 서버에서 pull:

```bash
cd /data_4tb/hifi-human-wgs-wdl-custom
git fetch --all && git reset --hard origin/main
chmod +x *.sh
```

> `git reset --hard`는 **로컬 변경사항을 모두 삭제**합니다. 중요한 수정이 있으면 먼저 백업.

### 특정 파일만 업데이트:

```bash
git fetch origin
git checkout origin/main -- create_batch_inputs.sh
git checkout origin/main -- batch_run_optimized.sh
```

## 자주 쓰는 명령어

```bash
git status                    # 변경사항 확인
git diff <file>               # 파일 변경 내용 확인
git log --oneline -10         # 최근 커밋 이력
git stash                     # 변경사항 임시 보관
git stash pop                 # 보관한 변경사항 복원
git branch -a                 # 전체 브랜치 목록
```

## 문제 해결

### .gitignore 규칙이 이미 추적된 파일에 적용 안 됨

```bash
git rm -r --cached .
git add .
git commit -m "Apply .gitignore"
```

### 실수로 대용량 파일 커밋 (아직 push 전)

```bash
git reset HEAD~1              # 마지막 커밋 취소 (파일은 유지)
```

### Conflict 발생 시

```bash
git reset --hard origin/main  # 원격 버전으로 완전히 덮어쓰기
```

### 원격 저장소 확인

```bash
git remote -v
git branch -r
```
