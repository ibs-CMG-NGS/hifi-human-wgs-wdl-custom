# Configuration Management Guide

## Input Configuration Files

### User-Specific Input Files (NOT tracked in git)
These files contain paths and settings specific to your local environment:
- `*.inputs.json` (e.g., `sample1.inputs.json`, `mysample.inputs.json`)
- These are in `.gitignore` to prevent tracking personal data paths

### Template Files (Tracked in git)
- `sample.inputs.json.example` - Template for singleton workflows
- `workflows/singleton.inputs.json` - Default template from the project
- `workflows/family.inputs.json` - Default template for family workflows

### How to Use
1. Copy the example file:
   ```bash
   cp sample.inputs.json.example my_sample.inputs.json
   ```

2. Edit `my_sample.inputs.json` with your specific:
   - Sample ID
   - Data file paths
   - Reference file locations
   - Backend configuration

3. Run the workflow:
   ```bash
   miniwdl run workflows/singleton.wdl -i my_sample.inputs.json
   ```

## Reference Files

### Tracked Templates (in git)
- `GRCh38.ref_map.v3p1p0.template.tsv`
- `GRCh38.tertiary_map.v3p1p0.template.tsv`

These are templates you can customize for your environment.

### Downloaded Resources (NOT in git)
- `hifi-wdl-resources/` - Large reference files downloaded separately
- These can be re-downloaded following the project documentation

## Environment Setup

### Tracked in git
- `environment.yml` - Conda environment specification
- `requirements.txt` - Python dependencies
- `setup_environment.sh` - Environment setup script
- `config/miniwdl.cfg` - Workflow engine configuration

### NOT tracked in git
- `.env` - Local environment variables
- `.venv/` - Python virtual environment directory

## Backend Configuration

The `backends/` directory contains backend-specific configurations:
- `backends/hpc/` - HPC cluster settings
- `backends/aws-healthomics/` - AWS configuration
- `backends/azure/` - Azure configuration
- `backends/gcp/` - GCP configuration

**Recommendation**: Create local copies of backend configs if you need to customize them:
```bash
cp backends/hpc/config.json backends/hpc/config.local.json
```

Then add `*.local.json` to your workflow to keep local customizations private.

## Best Practices

1. **Never commit**:
   - Personal data file paths
   - API keys or credentials
   - Large data files or reference genomes
   - Execution outputs and logs

2. **Always commit**:
   - Template/example configuration files
   - Documentation
   - Environment specifications
   - Workflow definition files (`.wdl`)

3. **For collaboration**:
   - Keep template files up-to-date
   - Document any required local customizations
   - Use relative paths when possible in templates
   - Document expected directory structure in README
