## 一套小流程，用于批处理Genbank文件，并提取蛋白序列

#### Step1 ----- 在NCBI Virus上Download所需accession，得到一个.acc文件

<img src="images/1.png" alt="NCBI Virus Accession_download" width="800">
<img src="images/2.png" alt="accession.acc" width="800">

#### Step2 ----- 将accession列表命名为accession_list.txt，编写批处理脚本，用于自动下载，详见download_gb_efetch.sh（与accession_list.txt放在同一目录下）
```bash
#!/usr/bin/env bash
# download_gb_efetch.sh
# 用法: ./download_gb_efetch.sh accession_list.txt out_dir [email] [api_key]

set -euo pipefail

if [[ $# -lt 2 ]]; then
  echo "用法: $0 accession_list.txt out_dir [email] [api_key]"
  exit 1
fi

ACC_FILE="$1"
OUT_DIR="$2"
EMAIL="${3:-}"
API_KEY="${4:-}"

mkdir -p "$OUT_DIR"

while IFS= read -r acc || [[ -n "$acc" ]]; do
  # 跳过空行和注释行
  [[ -z "$acc" || "$acc" =~ ^# ]] && continue
  out="$OUT_DIR/${acc}.gb"
  if [[ -s "$out" ]]; then
    echo "[跳过] $out 已存在且非空"
    continue
  fi
  url="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
  params="db=nuccore&id=${acc}&rettype=gb&retmode=text"
  [[ -n "$EMAIL" ]] && params="${params}&email=${EMAIL}"
  [[ -n "$API_KEY" ]] && params="${params}&api_key=${API_KEY}"

  echo "[下载] ${acc} -> ${out}"
  for attempt in {1..5}; do
    http_code=$(curl -sSL --retry 5 --retry-delay 2 --connect-timeout 10 \
      -w "%{http_code}" -o "$out".tmp "${url}?${params}" || true)
    if [[ "$http_code" == "200" ]] && grep -q "^LOCUS" "$out".tmp; then
      mv "$out".tmp "$out"
      break
    else
      echo "  尝试 ${attempt} 失败 (HTTP ${http_code})"
      sleep $((attempt*2))
    fi
    if [[ "$attempt" -eq 5 ]]; then
      echo "  错误: 放弃 ${acc}"; rm -f "$out".tmp
    fi
  done
done < "$ACC_FILE"

echo "下载完成"
```

#### Step3 ----- 在wsl中给脚本添加可执行权限
```bash
chmod +x download_gb_efetch.sh
```
