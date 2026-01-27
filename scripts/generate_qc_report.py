#!/usr/bin/env python3
"""
HiFi-human-WGS Pipeline QC Report Generator
여러 샘플의 분석 결과를 종합하여 HTML 리포트 생성
"""

import os
import json
import glob
from datetime import datetime
from pathlib import Path
import argparse

def parse_args():
    """명령줄 인자 파싱"""
    parser = argparse.ArgumentParser(
        description='Generate comprehensive QC report for HiFi WGS analysis'
    )
    parser.add_argument(
        '--batch-results',
        type=str,
        default='/data_4tb/hifi-human-wgs-wdl-custom/batch_results',
        help='Path to batch_results directory'
    )
    parser.add_argument(
        '--output',
        type=str,
        default='/data_4tb/hifi-human-wgs-wdl-custom/batch_results/QC_Report.html',
        help='Output HTML file path'
    )
    parser.add_argument(
        '--samples',
        nargs='+',
        help='Specific samples to include (default: all)'
    )
    return parser.parse_args()

def get_sample_names(batch_results_dir):
    """batch_results에서 샘플 이름 추출"""
    samples = []
    
    if not os.path.exists(batch_results_dir):
        print(f"Warning: Directory not found: {batch_results_dir}")
        return samples
    
    for item in os.listdir(batch_results_dir):
        item_path = os.path.join(batch_results_dir, item)
        if os.path.isdir(item_path) and item != 'logs':
            # _LAST 심볼릭 링크 또는 워크플로우 디렉토리가 있는지 확인
            last_link = os.path.join(item_path, '_LAST')
            # outputs.json이 있는 첫 번째 디렉토리 찾기
            workflow_dirs = [d for d in os.listdir(item_path) 
                           if os.path.isdir(os.path.join(item_path, d)) 
                           and d.startswith('202')]  # 날짜로 시작하는 디렉토리
            
            if os.path.exists(last_link) or workflow_dirs:
                samples.append(item)
    
    return sorted(samples)

def get_workflow_dir(sample_dir):
    """샘플 디렉토리에서 워크플로우 디렉토리 경로 반환 (_LAST 또는 최신 디렉토리)"""
    last_link = os.path.join(sample_dir, '_LAST')
    
    # _LAST 심볼릭 링크가 있으면 사용
    if os.path.exists(last_link):
        return last_link
    
    # 없으면 날짜로 시작하는 디렉토리 중 최신 것 찾기
    workflow_dirs = []
    for item in os.listdir(sample_dir):
        item_path = os.path.join(sample_dir, item)
        if os.path.isdir(item_path) and item.startswith('202'):  # 날짜 형식
            workflow_dirs.append((item, item_path))
    
    if workflow_dirs:
        # 디렉토리명 기준 정렬 (최신이 마지막)
        workflow_dirs.sort()
        return workflow_dirs[-1][1]
    
    return None

def parse_outputs_json(sample_dir):
    """outputs.json 파싱"""
    workflow_dir = get_workflow_dir(sample_dir)
    if not workflow_dir:
        return None
    
    outputs_file = os.path.join(workflow_dir, 'outputs.json')
    if not os.path.exists(outputs_file):
        return None
    
    try:
        with open(outputs_file, 'r') as f:
            return json.load(f)
    except Exception as e:
        print(f"Warning: Could not parse {outputs_file}: {e}")
        return None

def parse_bam_stats(sample_dir):
    """BAM statistics 파싱 - QC metrics 포함"""
    workflow_dir = get_workflow_dir(sample_dir)
    if not workflow_dir:
        return None
    
    stats_file = os.path.join(workflow_dir, 'out', 'bam_statistics', '*_stats.txt')
    stats_files = glob.glob(stats_file)
    
    if not stats_files:
        return None
    
    stats = {}
    try:
        with open(stats_files[0], 'r') as f:
            for line in f:
                if ':' in line:
                    key, value = line.strip().split(':', 1)
                    stats[key.strip()] = value.strip()
        
        # QC metrics 추출
        qc_metrics = {}
        
        # Read length (평균)
        if 'average length' in stats:
            try:
                qc_metrics['mean_read_length'] = float(stats['average length'])
            except:
                pass
        
        # Read quality (평균)
        if 'average quality' in stats:
            try:
                qc_metrics['mean_read_quality'] = float(stats['average quality'])
            except:
                pass
        
        # Total reads
        if 'sequences' in stats or 'total length' in stats:
            try:
                qc_metrics['total_reads'] = int(stats.get('sequences', '0').replace(',', ''))
            except:
                pass
        
        # Mapping rate (pbmm2 통계에서)
        if 'mapped' in stats:
            try:
                mapped = int(stats['mapped'].replace(',', ''))
                total = qc_metrics.get('total_reads', 0)
                if total > 0:
                    qc_metrics['mapping_rate'] = (mapped / total) * 100
            except:
                pass
        
        stats['qc_metrics'] = qc_metrics
        
    except Exception as e:
        print(f"Warning: Could not parse BAM stats: {e}")
        return None
    
    return stats

def parse_mosdepth_summary(sample_dir):
    """mosdepth summary 파싱"""
    workflow_dir = get_workflow_dir(sample_dir)
    if not workflow_dir:
        return {'mean_coverage': 0}
    
    summary_file = os.path.join(workflow_dir, 'out', 'mosdepth_summary', '*.mosdepth.summary.txt')
    summary_files = glob.glob(summary_file)
    
    if not summary_files:
        return {'mean_coverage': 0}
    
    data = {}
    try:
        with open(summary_files[0], 'r') as f:
            lines = f.readlines()
            # 헤더 스킵
            for line in lines[1:]:  # 첫 번째 줄은 헤더
                if line.startswith('total') or line.startswith('chr'):
                    parts = line.strip().split('\t')
                    if len(parts) >= 4:
                        region = parts[0]
                        try:
                            mean_depth = float(parts[3])
                            if region == 'total':
                                data['mean_coverage'] = mean_depth
                                break
                        except ValueError:
                            continue
    except Exception as e:
        print(f"Warning: Could not parse mosdepth summary: {e}")
        return {'mean_coverage': 0}
    
    # mean_coverage가 없으면 기본값 설정
    if 'mean_coverage' not in data:
        data['mean_coverage'] = 0
    
    return data

def parse_workflow_log(sample_dir):
    """workflow.log에서 실행 시간 정보 파싱"""
    workflow_dir = get_workflow_dir(sample_dir)
    if not workflow_dir:
        return None
    
    log_file = os.path.join(workflow_dir, 'workflow.log')
    if not os.path.exists(log_file):
        return None
    
    data = {
        'start_time': None,
        'end_time': None,
        'duration_seconds': 0,
        'duration_formatted': 'Unknown'
    }
    
    try:
        with open(log_file, 'r') as f:
            lines = f.readlines()
            
            # 첫 줄에서 시작 시간 추출
            if lines:
                first_line = lines[0]
                # Format: "2026-01-20 10:17:04.123 ..."
                import re
                match = re.search(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', first_line)
                if match:
                    data['start_time'] = match.group(1)
            
            # 마지막 줄에서 종료 시간 추출
            if len(lines) > 1:
                last_line = lines[-1]
                match = re.search(r'(\d{4}-\d{2}-\d{2} \d{2}:\d{2}:\d{2})', last_line)
                if match:
                    data['end_time'] = match.group(1)
            
            # 실행 시간 계산
            if data['start_time'] and data['end_time']:
                from datetime import datetime
                start = datetime.strptime(data['start_time'], '%Y-%m-%d %H:%M:%S')
                end = datetime.strptime(data['end_time'], '%Y-%m-%d %H:%M:%S')
                duration = end - start
                data['duration_seconds'] = int(duration.total_seconds())
                
                # 포맷팅
                hours = data['duration_seconds'] // 3600
                minutes = (data['duration_seconds'] % 3600) // 60
                data['duration_formatted'] = f"{hours}h {minutes}m"
    
    except Exception as e:
        print(f"Warning: Could not parse workflow log: {e}")
        return None
    
    return data

def parse_small_variant_stats(sample_dir):
    """Small variant statistics 파싱"""
    workflow_dir = get_workflow_dir(sample_dir)
    if not workflow_dir:
        return None
    
    stats_file = os.path.join(workflow_dir, 'out', 'small_variant_stats', '*.stats.txt')
    stats_files = glob.glob(stats_file)
    
    if not stats_files:
        return None
    
    data = {}
    try:
        with open(stats_files[0], 'r') as f:
            content = f.read()
            # SNP 개수 추출
            if 'number of SNPs:' in content:
                for line in content.split('\n'):
                    if 'number of SNPs:' in line:
                        data['snps'] = int(line.split(':')[1].strip())
            # Indel 개수 추출
            if 'number of indels:' in content:
                for line in content.split('\n'):
                    if 'number of indels:' in line:
                        data['indels'] = int(line.split(':')[1].strip())
    except Exception as e:
        print(f"Warning: Could not parse variant stats: {e}")
        return None
    
    return data

def get_file_size(filepath):
    """파일 크기를 읽기 쉬운 형식으로 변환"""
    if not os.path.exists(filepath):
        return "N/A"
    
    size = os.path.getsize(filepath)
    if size > 1e12:
        return f"{size/1e12:.2f} TB"
    elif size > 1e9:
        return f"{size/1e9:.2f} GB"
    elif size > 1e6:
        return f"{size/1e6:.2f} MB"
    else:
        return f"{size/1e3:.2f} KB"

def count_vcf_variants(vcf_file):
    """VCF 파일의 variant 개수 카운트"""
    if not os.path.exists(vcf_file):
        return 0
    
    import subprocess
    try:
        # bcftools가 있으면 사용
        result = subprocess.run(
            ['bcftools', 'view', '-H', vcf_file, '|', 'wc', '-l'],
            shell=True,
            capture_output=True,
            text=True
        )
        return int(result.stdout.strip())
    except:
        # bcftools가 없으면 직접 카운트
        count = 0
        try:
            import gzip
            with gzip.open(vcf_file, 'rt') as f:
                for line in f:
                    if not line.startswith('#'):
                        count += 1
        except:
            pass
        return count

def collect_sample_data(batch_results_dir, sample):
    """샘플별 데이터 수집"""
    sample_dir = os.path.join(batch_results_dir, sample)
    
    data = {
        'sample_name': sample,
        'outputs': parse_outputs_json(sample_dir),
        'bam_stats': parse_bam_stats(sample_dir),
        'coverage': parse_mosdepth_summary(sample_dir),
        'variant_stats': parse_small_variant_stats(sample_dir),
        'workflow_log': parse_workflow_log(sample_dir)
    }
    
    # 주요 파일 크기
    workflow_dir = get_workflow_dir(sample_dir)
    if workflow_dir:
        out_dir = os.path.join(workflow_dir, 'out')
    else:
        out_dir = ''
    
    data['file_sizes'] = {
        'aligned_bam': get_file_size(
            glob.glob(os.path.join(out_dir, 'merged_haplotagged_bam', '*.bam'))[0]
            if out_dir and glob.glob(os.path.join(out_dir, 'merged_haplotagged_bam', '*.bam')) else ''
        ),
        'small_variant_vcf': get_file_size(
            glob.glob(os.path.join(out_dir, 'phased_small_variant_vcf', '*.vcf.gz'))[0]
            if out_dir and glob.glob(os.path.join(out_dir, 'phased_small_variant_vcf', '*.vcf.gz')) else ''
        ),
        'sv_vcf': get_file_size(
            glob.glob(os.path.join(out_dir, 'phased_sv_vcf', '*.vcf.gz'))[0]
            if out_dir and glob.glob(os.path.join(out_dir, 'phased_sv_vcf', '*.vcf.gz')) else ''
        )
    }
    
    # Variant 개수
    small_vcf = glob.glob(os.path.join(out_dir, 'phased_small_variant_vcf', '*.vcf.gz')) if out_dir else []
    sv_vcf = glob.glob(os.path.join(out_dir, 'phased_sv_vcf', '*.vcf.gz')) if out_dir else []
    
    data['variant_counts'] = {
        'small_variants': count_vcf_variants(small_vcf[0]) if small_vcf else 0,
        'structural_variants': count_vcf_variants(sv_vcf[0]) if sv_vcf else 0
    }
    
    # 로그 파일에서 완료 시간 추출
    log_file = os.path.join(batch_results_dir, 'logs', f'{sample}.log')
    if os.path.exists(log_file):
        try:
            mtime = os.path.getmtime(log_file)
            data['completion_time'] = datetime.fromtimestamp(mtime).strftime('%Y-%m-%d %H:%M:%S')
        except:
            data['completion_time'] = 'Unknown'
    else:
        data['completion_time'] = 'Unknown'
    
    return data

def generate_html_report(samples_data, batch_results_dir):
    """HTML 리포트 생성"""
    
    total_samples = len(samples_data)
    
    # 전체 통계 계산
    total_coverage = sum(
        (s.get('coverage') or {}).get('mean_coverage', 0)
        for s in samples_data 
        if s.get('coverage') and isinstance(s.get('coverage'), dict)
    )
    avg_coverage = total_coverage / total_samples if total_samples > 0 else 0
    
    total_variants = sum(
        s['variant_counts'].get('small_variants', 0) 
        for s in samples_data
    )
    total_svs = sum(
        s['variant_counts'].get('structural_variants', 0) 
        for s in samples_data
    )
    
    html = f"""
<!DOCTYPE html>
<html lang="ko">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>HiFi WGS Pipeline QC Report</title>
    <script src="https://cdn.jsdelivr.net/npm/chart.js@4.4.0/dist/chart.umd.min.js"></script>
    <style>
        * {{
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }}
        
        body {{
            font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            padding: 20px;
            color: #333;
        }}
        
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            background: white;
            border-radius: 15px;
            box-shadow: 0 10px 40px rgba(0,0,0,0.2);
            overflow: hidden;
        }}
        
        .header {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 40px;
            text-align: center;
        }}
        
        .header h1 {{
            font-size: 2.5em;
            margin-bottom: 10px;
        }}
        
        .header p {{
            font-size: 1.1em;
            opacity: 0.9;
        }}
        
        .content {{
            padding: 40px;
        }}
        
        .section {{
            margin-bottom: 40px;
        }}
        
        .section-title {{
            font-size: 1.8em;
            color: #667eea;
            margin-bottom: 20px;
            padding-bottom: 10px;
            border-bottom: 3px solid #667eea;
        }}
        
        .summary-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(250px, 1fr));
            gap: 20px;
            margin-bottom: 30px;
        }}
        
        .summary-card {{
            background: linear-gradient(135deg, #f5f7fa 0%, #c3cfe2 100%);
            padding: 25px;
            border-radius: 10px;
            box-shadow: 0 4px 6px rgba(0,0,0,0.1);
            transition: transform 0.3s ease;
        }}
        
        .summary-card:hover {{
            transform: translateY(-5px);
        }}
        
        .summary-card h3 {{
            color: #667eea;
            font-size: 0.9em;
            margin-bottom: 10px;
            text-transform: uppercase;
        }}
        
        .summary-card .value {{
            font-size: 2em;
            font-weight: bold;
            color: #333;
        }}
        
        .summary-card .sub-value {{
            font-size: 0.9em;
            color: #666;
            margin-top: 5px;
        }}
        
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-top: 20px;
            background: white;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
            border-radius: 8px;
            overflow: hidden;
        }}
        
        th {{
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 15px;
            text-align: left;
            font-weight: 600;
            position: sticky;
            top: 0;
        }}
        
        td {{
            padding: 12px 15px;
            border-bottom: 1px solid #eee;
        }}
        
        tr:hover {{
            background-color: #f8f9ff;
        }}
        
        .metric-good {{
            color: #10b981;
            font-weight: bold;
        }}
        
        .metric-warning {{
            color: #f59e0b;
            font-weight: bold;
        }}
        
        .metric-bad {{
            color: #ef4444;
            font-weight: bold;
        }}
        
        .badge {{
            display: inline-block;
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 0.85em;
            font-weight: 600;
        }}
        
        .badge-success {{
            background: #d1fae5;
            color: #065f46;
        }}
        
        .badge-warning {{
            background: #fef3c7;
            color: #92400e;
        }}
        
        .badge-info {{
            background: #dbeafe;
            color: #1e40af;
        }}
        
        .progress-bar {{
            height: 25px;
            background: #e5e7eb;
            border-radius: 12px;
            overflow: hidden;
            margin: 5px 0;
        }}
        
        .progress-fill {{
            height: 100%;
            background: linear-gradient(90deg, #667eea 0%, #764ba2 100%);
            display: flex;
            align-items: center;
            justify-content: center;
            color: white;
            font-size: 0.85em;
            font-weight: bold;
            transition: width 0.3s ease;
        }}
        
        .footer {{
            background: #f9fafb;
            padding: 20px;
            text-align: center;
            color: #666;
            font-size: 0.9em;
        }}
        
        .sample-header {{
            background: #f0f4ff;
            padding: 15px;
            margin: 20px 0;
            border-radius: 8px;
            border-left: 4px solid #667eea;
        }}
        
        .sample-header h3 {{
            color: #667eea;
            margin-bottom: 5px;
        }}
        
        .chart-container {{
            background: white;
            padding: 30px;
            border-radius: 10px;
            margin: 20px 0;
            box-shadow: 0 2px 10px rgba(0,0,0,0.1);
        }}
        
        .chart-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(400px, 1fr));
            gap: 20px;
            margin: 20px 0;
        }}
        
        canvas {{
            max-height: 400px;
        }}
        
        .stats-box {{
            background: #f8f9ff;
            padding: 15px;
            border-radius: 8px;
            margin: 10px 0;
        }}
        
        .stats-box h4 {{
            color: #667eea;
            margin-bottom: 10px;
            font-size: 1.1em;
        }}
        
        .stats-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
            gap: 15px;
        }}
        
        .stat-item {{
            background: white;
            padding: 12px;
            border-radius: 6px;
            border-left: 3px solid #667eea;
        }}
        
        .stat-label {{
            font-size: 0.85em;
            color: #666;
            margin-bottom: 5px;
        }}
        
        .stat-value {{
            font-size: 1.3em;
            font-weight: bold;
            color: #333;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>🧬 HiFi WGS Pipeline QC Report</h1>
            <p>Comprehensive Quality Control & Analysis Summary</p>
            <p style="font-size: 0.9em; margin-top: 10px;">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
        
        <div class="content">
            <!-- 전체 요약 -->
            <div class="section">
                <h2 class="section-title">📊 Overall Summary</h2>
                <div class="summary-grid">
                    <div class="summary-card">
                        <h3>Total Samples</h3>
                        <div class="value">{total_samples}</div>
                        <div class="sub-value">Successfully analyzed</div>
                    </div>
                    <div class="summary-card">
                        <h3>Average Coverage</h3>
                        <div class="value">{avg_coverage:.1f}×</div>
                        <div class="sub-value">Genome-wide mean depth</div>
                    </div>
                    <div class="summary-card">
                        <h3>Total Variants</h3>
                        <div class="value">{total_variants:,}</div>
                        <div class="sub-value">Small variants (SNPs + Indels)</div>
                    </div>
                    <div class="summary-card">
                        <h3>Structural Variants</h3>
                        <div class="value">{total_svs:,}</div>
                        <div class="sub-value">SVs detected across samples</div>
                    </div>
                </div>
            </div>
            
            <!-- 샘플 목록 -->
            <div class="section">
                <h2 class="section-title">📋 Sample List</h2>
                <div class="summary-grid">
"""
    
    for data in samples_data:
        sample = data['sample_name']
        completion = data.get('completion_time', 'Unknown')
        html += f"""
                    <div class="summary-card">
                        <h3>{sample}</h3>
                        <div class="sub-value">Completed: {completion}</div>
                    </div>
"""
    
    html += """
                </div>
            </div>
            
            <!-- 시각화 차트 -->
            <div class="section">
                <h2 class="section-title">📊 Visual Analytics</h2>
                <div class="chart-grid">
                    <div class="chart-container">
                        <canvas id="coverageChart"></canvas>
                    </div>
                    <div class="chart-container">
                        <canvas id="variantChart"></canvas>
                    </div>
                </div>
            </div>
            
            <!-- Coverage 통계 -->
            <div class="section">
                <h2 class="section-title">📈 Coverage Statistics</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Mean Coverage</th>
                            <th>Coverage Quality</th>
                            <th>Aligned BAM Size</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    for data in samples_data:
        sample = data['sample_name']
        coverage_data = data.get('coverage') or {}
        coverage = coverage_data.get('mean_coverage', 0) if isinstance(coverage_data, dict) else 0
        bam_size = data['file_sizes'].get('aligned_bam', 'N/A')
        
        # Coverage 평가
        if coverage >= 30:
            badge = '<span class="badge badge-success">Excellent</span>'
            color_class = 'metric-good'
        elif coverage >= 20:
            badge = '<span class="badge badge-warning">Good</span>'
            color_class = 'metric-warning'
        else:
            badge = '<span class="badge badge-warning">Low</span>'
            color_class = 'metric-warning'
        
        html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td class="{color_class}">{coverage:.1f}×</td>
                            <td>{badge}</td>
                            <td>{bam_size}</td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- Variant Calling 결과 -->
            <div class="section">
                <h2 class="section-title">🧪 Variant Calling Results</h2>
                <table>
                    <thead>
                        <tr>
                            <th>Sample</th>
                            <th>Small Variants</th>
                            <th>Structural Variants</th>
                            <th>Small VCF Size</th>
                            <th>SV VCF Size</th>
                        </tr>
                    </thead>
                    <tbody>
"""
    
    for data in samples_data:
        sample = data['sample_name']
        small_vars = data['variant_counts'].get('small_variants', 0)
        svs = data['variant_counts'].get('structural_variants', 0)
        small_vcf_size = data['file_sizes'].get('small_variant_vcf', 'N/A')
        sv_vcf_size = data['file_sizes'].get('sv_vcf', 'N/A')
        
        html += f"""
                        <tr>
                            <td><strong>{sample}</strong></td>
                            <td class="metric-good">{small_vars:,}</td>
                            <td class="metric-info">{svs:,}</td>
                            <td>{small_vcf_size}</td>
                            <td>{sv_vcf_size}</td>
                        </tr>
"""
    
    html += """
                    </tbody>
                </table>
            </div>
            
            <!-- 상세 샘플 정보 -->
            <div class="section">
                <h2 class="section-title">🔍 Detailed Sample Information</h2>
"""
    
    for data in samples_data:
        sample = data['sample_name']
        outputs = data.get('outputs', {})
        
        html += f"""
                <div class="sample-header">
                    <h3>{sample}</h3>
                    <p>Completion time: {data.get('completion_time', 'Unknown')}</p>
                </div>
                
                <table>
                    <thead>
                        <tr>
                            <th>Output Type</th>
                            <th>File Path</th>
                            <th>Status</th>
                        </tr>
                    </thead>
                    <tbody>
"""
        
        # 주요 출력 파일 정보
        key_outputs = [
            ('humanwgs_singleton.phased_small_variant_vcf', 'Phased Small Variant VCF'),
            ('humanwgs_singleton.phased_sv_vcf', 'Phased SV VCF'),
            ('humanwgs_singleton.merged_haplotagged_bam', 'Haplotagged BAM'),
            ('humanwgs_singleton.pharmcat_report_html', 'PharmCAT Report'),
            ('humanwgs_singleton.phase_stats', 'Phasing Statistics'),
            ('humanwgs_singleton.mosdepth_summary', 'Coverage Summary'),
            ('humanwgs_singleton.bam_statistics', 'BAM Statistics')
        ]
        
        if outputs:
            for key, label in key_outputs:
                if key in outputs:
                    value = outputs[key]
                    if isinstance(value, list) and value:
                        value = value[0]
                    
                    file_exists = os.path.exists(value) if isinstance(value, str) else False
                    status = '<span class="badge badge-success">✓ Available</span>' if file_exists else '<span class="badge badge-warning">⚠ Not Found</span>'
                    filename = os.path.basename(value) if isinstance(value, str) else 'N/A'
                    
                    html += f"""
                        <tr>
                            <td><strong>{label}</strong></td>
                            <td style="font-family: monospace; font-size: 0.9em;">{filename}</td>
                            <td>{status}</td>
                        </tr>
"""
        
        html += """
                    </tbody>
                </table>
"""
    
    html += f"""
            </div>
        </div>
        
        <div class="footer">
            <p>Generated by HiFi WGS Pipeline QC Report Generator</p>
            <p>Pipeline: PacBio HiFi-human-WGS-WDL v3.1.0</p>
            <p>Report Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}</p>
        </div>
    </div>
    
    <script>
        // Coverage 차트 데이터 준비
        const coverageData = {{
            labels: {json.dumps([d['sample_name'] for d in samples_data])},
            datasets: [{{
                label: 'Mean Coverage (×)',
                data: {json.dumps([(d.get('coverage') or {{}}).get('mean_coverage', 0) for d in samples_data])},
                backgroundColor: 'rgba(102, 126, 234, 0.6)',
                borderColor: 'rgba(102, 126, 234, 1)',
                borderWidth: 2
            }}]
        }};
        
        // Variant 차트 데이터 준비
        const variantData = {{
            labels: {json.dumps([d['sample_name'] for d in samples_data])},
            datasets: [
                {{
                    label: 'Small Variants',
                    data: {json.dumps([d['variant_counts'].get('small_variants', 0) for d in samples_data])},
                    backgroundColor: 'rgba(16, 185, 129, 0.6)',
                    borderColor: 'rgba(16, 185, 129, 1)',
                    borderWidth: 2
                }},
                {{
                    label: 'Structural Variants',
                    data: {json.dumps([d['variant_counts'].get('structural_variants', 0) for d in samples_data])},
                    backgroundColor: 'rgba(245, 158, 11, 0.6)',
                    borderColor: 'rgba(245, 158, 11, 1)',
                    borderWidth: 2
                }}
            ]
        }};
        
        // Coverage 차트 생성
        const coverageCtx = document.getElementById('coverageChart').getContext('2d');
        new Chart(coverageCtx, {{
            type: 'bar',
            data: coverageData,
            options: {{
                responsive: true,
                maintainAspectRatio: true,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Coverage Distribution Across Samples',
                        font: {{ size: 16, weight: 'bold' }}
                    }},
                    legend: {{ display: true }}
                }},
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{ display: true, text: 'Coverage (×)' }}
                    }}
                }}
            }}
        }});
        
        // Variant 차트 생성
        const variantCtx = document.getElementById('variantChart').getContext('2d');
        new Chart(variantCtx, {{
            type: 'bar',
            data: variantData,
            options: {{
                responsive: true,
                maintainAspectRatio: true,
                plugins: {{
                    title: {{
                        display: true,
                        text: 'Variant Counts Across Samples',
                        font: {{ size: 16, weight: 'bold' }}
                    }},
                    legend: {{ display: true, position: 'top' }}
                }},
                scales: {{
                    y: {{
                        beginAtZero: true,
                        title: {{ display: true, text: 'Number of Variants' }}
                    }}
                }}
            }}
        }});
    </script>
</body>
</html>
"""
    
    return html

def main():
    """메인 실행 함수"""
    args = parse_args()
    
    print("🧬 HiFi WGS Pipeline QC Report Generator")
    print("=" * 60)
    print(f"Batch results directory: {args.batch_results}")
    print(f"Output file: {args.output}")
    
    # 샘플 이름 수집
    if args.samples:
        samples = args.samples
        print(f"✓ Using specified {len(samples)} sample(s): {', '.join(samples)}")
    else:
        samples = get_sample_names(args.batch_results)
        print(f"✓ Found {len(samples)} sample(s): {', '.join(samples)}")
    
    if not samples:
        print("❌ No samples found!")
        return 1
    
    # 각 샘플 데이터 수집
    print("\n📊 Collecting sample data...")
    samples_data = []
    for sample in samples:
        print(f"  - Processing {sample}...")
        data = collect_sample_data(args.batch_results, sample)
        samples_data.append(data)
    
    # HTML 리포트 생성
    print("\n📝 Generating HTML report...")
    html_content = generate_html_report(samples_data, args.batch_results)
    
    # 파일 저장
    output_dir = os.path.dirname(args.output)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html_content)
    
    print(f"✅ Report saved: {args.output}")
    print(f"\n🌐 Open the report in your browser:")
    print(f"   file://{os.path.abspath(args.output)}")
    
    # 터미널 요약
    print("\n" + "=" * 60)
    print("📊 Quick Summary")
    print("=" * 60)
    
    for data in samples_data:
        sample = data['sample_name']
        coverage_data = data.get('coverage') or {}
        coverage = coverage_data.get('mean_coverage', 0) if isinstance(coverage_data, dict) else 0
        small_vars = data['variant_counts'].get('small_variants', 0)
        svs = data['variant_counts'].get('structural_variants', 0)
        
        print(f"\n{sample}:")
        print(f"  Coverage: {coverage:.1f}×")
        print(f"  Small variants: {small_vars:,}")
        print(f"  Structural variants: {svs:,}")
    
    return 0

if __name__ == '__main__':
    exit(main())
