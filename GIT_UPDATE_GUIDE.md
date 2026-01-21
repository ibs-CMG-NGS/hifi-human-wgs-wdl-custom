# Git으로 Batch Processing 스크립트 업데이트하기

## 📤 Windows/VS Code에서 Git Push

### 1단계: Git 상태 확인
```powershell
cd \\wsl.localhost\Ubuntu\home\ygkim\ngs_pipeline\HiFi-human-WGS-WDL
git status
```

### 2단계: 새로 생성된 파일 추가
```powershell
git add create_batch_inputs.sh
git add batch_run_optimized.sh
git add monitor_batch.sh
git add collect_results.sh
git add samples.csv
git add BATCH_PROCESSING_GUIDE.md
git add BATCH_QUICK_START.md
git add GIT_UPDATE_GUIDE.md
```

또는 한 번에:
```powershell
git add *.sh *.csv *.md
```

### 3단계: 커밋
```powershell
git commit -m "Add batch processing scripts and guides"
```

### 4단계: Push
```powershell
git push origin main
```

또는 브랜치가 다르면:
```powershell
git push origin master
```

---

## 📥 서버에서 Git Pull (강제 덮어쓰기)

서버에 SSH 접속 후:

### 방법 1: 로컬 변경사항 완전히 무시하고 덮어쓰기 (권장)

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# 현재 상태 확인
git status

# 로컬 변경사항 모두 삭제 (주의!)
git reset --hard HEAD

# 원격 저장소 최신 상태 가져오기
git fetch origin

# 강제로 원격 브랜치로 덮어쓰기
git reset --hard origin/main
# 또는 master 브랜치인 경우:
# git reset --hard origin/master

# Pull 실행
git pull origin main
```

### 방법 2: 특정 파일만 원격 저장소 버전으로 교체

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom

# 원격 저장소 최신 정보 가져오기
git fetch origin

# 특정 파일만 원격 저장소 버전으로 교체
git checkout origin/main -- create_batch_inputs.sh
git checkout origin/main -- batch_run_optimized.sh
git checkout origin/main -- monitor_batch.sh
git checkout origin/main -- collect_results.sh
git checkout origin/main -- samples.csv
git checkout origin/main -- BATCH_PROCESSING_GUIDE.md
git checkout origin/main -- BATCH_QUICK_START.md
```

### 방법 3: 한 줄 명령어 (가장 간단)

```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom
git fetch --all && git reset --hard origin/main && git pull origin main
```

---

## ✅ 업데이트 후 확인

### 1. 파일 존재 확인
```bash
ls -lh *.sh *.csv *BATCH*.md
```

### 2. 실행 권한 부여
```bash
chmod +x create_batch_inputs.sh
chmod +x batch_run_optimized.sh
chmod +x monitor_batch.sh
chmod +x collect_results.sh
```

### 3. 파일 내용 확인
```bash
head -5 batch_run_optimized.sh
```

---

## 🔄 전체 업데이트 프로세스 요약

### Windows에서:
```powershell
cd \\wsl.localhost\Ubuntu\home\ygkim\ngs_pipeline\HiFi-human-WGS-WDL
git add *.sh *.csv *BATCH*.md GIT_UPDATE_GUIDE.md
git commit -m "Add batch processing scripts"
git push origin main
```

### 서버에서:
```bash
cd ~/ngs-pipeline/hifi-human-wgs-wdl-custom
git fetch --all && git reset --hard origin/main && git pull origin main
chmod +x *.sh
```

---

## ⚠️ 주의사항

1. **`git reset --hard`는 로컬 변경사항을 모두 삭제합니다!**
   - 중요한 수정사항이 있다면 백업하세요
   - 또는 방법 2를 사용하세요

2. **브랜치 이름 확인**
   - `main` 또는 `master` 브랜치 이름 확인 필요
   ```bash
   git branch -a
   ```

3. **권한 문제**
   - Pull 후 스크립트 실행 권한 재설정 필요
   ```bash
   chmod +x *.sh
   ```

4. **samples.csv는 서버에서 수정 필요**
   - Git으로 받은 후 실제 샘플 정보로 편집
   ```bash
   vim samples.csv
   ```

---

## 🚨 문제 해결

### "Permission denied" 에러
```bash
chmod +x *.sh
```

### "conflict" 에러
```bash
git reset --hard origin/main
git pull origin main
```

### 브랜치 이름 모를 때
```bash
git branch -r
# origin/main 또는 origin/master 확인
```

### 원격 저장소 URL 확인
```bash
git remote -v
```

---

이제 Windows에서 Git push만 하면, 서버에서 한 줄 명령어로 최신 스크립트를 받을 수 있습니다!
