#!/usr/bin/env bash
set -eo pipefail

# WDL 이 참조하는 docker 이미지 전부를 manifest 로 뽑는다.
# 출력은 ./image_manifest.txt, 소비자는 scripts/populate_miniwdl_singularity_cache.sh.
#
# 2026-11-26: **태그 기반 이미지를 놓치고 있었다.** 원래는 `grep '@sha'` 만
# 해서 quay.io/pacbio 의 digest 고정 이미지 16개만 잡혔고, biocontainers·
# ontresearch 처럼 `name:tag` 로 적힌 16개는 manifest 에 안 들어갔다.
# (그중 일부가 과거에 손으로 덧붙여져 있어서 결손이 더 안 보였다.)
# manifest 가 불완전하면 캐시를 미리 채워도 첫 실행이 빠진 이미지를 받으러
# 나간다 — 그게 오프라인/혼잡한 노드에서 실행을 세운다.

OUT=./image_manifest.txt

# (1) digest 고정 이미지. container_registry 플레이스홀더를 기본값으로 치환한다.
grep '@sha' -h -r workflows/ \
| tr --squeeze-repeats ' ' \
| cut --fields=3 --delimiter=' ' \
| sed 's!~{runtime_attributes.container_registry}!quay.io/pacbio!;s/"//g;' \
| sort --unique \
> "${OUT}"

# (2) 태그 기반 이미지. `local/` 은 레지스트리에 없는 로컬 빌드라 제외한다
#     (local/kineticstools:0.8.0 — 쓰려면 직접 빌드해야 한다).
grep -rh 'docker\|container' workflows/ --include=*.wdl \
| grep -oE '"[a-z0-9._-]+(\.[a-z]+)?/[A-Za-z0-9._/-]+:[A-Za-z0-9._-]+"' \
| tr -d '"' \
| grep -v '^local/' \
| sort --unique \
>> "${OUT}"

# (3) 워크플로가 버전 변수로 조립하는 이미지. backend != AWS-HealthOmics 이면
#     deepvariant 레지스트리는 google 로 고정되고, gpu=true 작업은 같은 태그에
#     "-gpu" 가 붙는다 (workflows/.../deepvariant.wdl 의 gpu 분기).
deepvariant_version=1.9.0
echo "google/deepvariant:${deepvariant_version}" >> "${OUT}"
echo "google/deepvariant:${deepvariant_version}-gpu" >> "${OUT}"

# run_pgx=true 일 때만 쓴다. 안 쓰면 받을 필요 없다.
pharmcat_version=2.15.4
echo "pgkb/pharmcat:${pharmcat_version}" >> "${OUT}"

sort -u -o "${OUT}" "${OUT}"
echo "$(wc -l < "${OUT}") images -> ${OUT}" >&2
