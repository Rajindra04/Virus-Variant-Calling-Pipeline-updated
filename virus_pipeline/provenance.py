"""
Provenance tracking for the Virus Variant Calling Pipeline.

Records every tool, version, parameters, and outcome for each pipeline step.
Outputs a JSON log and human-readable report for reproducibility and auditing.
"""

import json
import os
import subprocess
import logging
import datetime

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')


class ProvenanceTracker:
    """Tracks all pipeline steps, tools, parameters, and outcomes."""

    def __init__(self, output_dir):
        self.output_dir = output_dir
        self.start_time = datetime.datetime.now().isoformat()
        self.steps = []
        self.tool_versions = {}
        self.pipeline_args = {}
        self.config = {}

    def set_pipeline_args(self, args_dict):
        """Record the top-level pipeline arguments."""
        self.pipeline_args = {k: str(v) for k, v in args_dict.items()}

    def set_config(self, config_dict):
        """Record the full virus config used for this run."""
        self.config = config_dict

    def detect_tool_version(self, tool_name, version_cmd):
        """Detect and record the version of a tool."""
        try:
            result = subprocess.run(
                version_cmd, shell=True,
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=30
            )
            output = (result.stdout.decode('utf-8').strip() or
                      result.stderr.decode('utf-8').strip())
            # Take just the first line for cleanliness
            version_str = output.split('\n')[0].strip()
            self.tool_versions[tool_name] = version_str
        except Exception as e:
            self.tool_versions[tool_name] = f"version detection failed: {e}"

    def detect_all_tool_versions(self):
        """Detect versions for all tools used in the pipeline."""
        version_commands = {
            'fastp': 'fastp --version 2>&1',
            'fastqc': 'fastqc --version 2>&1',
            'bwa-mem2': 'bwa-mem2 version 2>&1',
            'samtools': 'samtools --version 2>&1 | head -1',
            'gatk': 'gatk --version 2>&1 | head -2 | tail -1',
            'ivar': 'ivar version 2>&1 | head -1',
            'snpEff': 'snpEff -version 2>&1 | head -1',
            'SnpSift': 'SnpSift -version 2>&1 | head -1',
        }
        for tool, cmd in version_commands.items():
            self.detect_tool_version(tool, cmd)

    def record_step(self, step_name, tool, parameters, status="completed",
                    sample=None, notes=None):
        """Record a single pipeline step.

        Args:
            step_name: Human-readable step name (e.g., "Read trimming")
            tool: Tool used (e.g., "fastp")
            parameters: Dict of parameters used
            status: "completed", "skipped", or "failed"
            sample: Sample name if per-sample step, None if global
            notes: Optional notes (e.g., why a step was skipped)
        """
        entry = {
            'timestamp': datetime.datetime.now().isoformat(),
            'step': step_name,
            'tool': tool,
            'parameters': parameters,
            'status': status,
            'sample': sample,
        }
        if notes:
            entry['notes'] = notes
        self.steps.append(entry)

        # Log it
        status_icon = {'completed': 'RAN', 'skipped': 'SKIPPED', 'failed': 'FAILED'}
        icon = status_icon.get(status, status.upper())
        sample_str = f" [{sample}]" if sample else ""
        logging.info(f"[PROVENANCE] {icon}: {step_name}{sample_str} ({tool})")
        if notes:
            logging.info(f"[PROVENANCE]   Note: {notes}")

    def to_dict(self):
        """Return the full provenance record as a dictionary."""
        result = {
            'pipeline': 'Virus-Variant-Calling-Pipeline',
            'run_start': self.start_time,
            'run_end': datetime.datetime.now().isoformat(),
            'pipeline_arguments': self.pipeline_args,
            'tool_versions': self.tool_versions,
            'steps': self.steps,
        }
        if self.config:
            result['virus_config'] = self.config
        return result

    def write_json(self, filename='provenance.json'):
        """Write provenance log as JSON."""
        path = os.path.join(self.output_dir, filename)
        with open(path, 'w') as f:
            json.dump(self.to_dict(), f, indent=2)
        logging.info(f"Provenance JSON written to {path}")
        return path

    def write_report(self, filename='provenance_report.txt'):
        """Write a human-readable provenance report."""
        path = os.path.join(self.output_dir, filename)
        data = self.to_dict()

        with open(path, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("VIRUS VARIANT CALLING PIPELINE — PROVENANCE REPORT\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"Run started:  {data['run_start']}\n")
            f.write(f"Run finished: {data['run_end']}\n\n")

            # Pipeline arguments
            f.write("-" * 80 + "\n")
            f.write("PIPELINE ARGUMENTS\n")
            f.write("-" * 80 + "\n")
            for k, v in data['pipeline_arguments'].items():
                f.write(f"  {k}: {v}\n")
            f.write("\n")

            # Virus config
            if data.get('virus_config'):
                vc = data['virus_config']
                f.write("-" * 80 + "\n")
                f.write("VIRUS CONFIGURATION\n")
                f.write("-" * 80 + "\n")
                f.write(f"  Virus: {vc.get('virus_name', 'N/A')}\n")
                f.write(f"  Genome type: {vc.get('genome_type', 'N/A')}\n")
                f.write(f"  Ploidy: {vc.get('ploidy', 'N/A')}\n")
                f.write(f"  Expected genome size: {vc.get('expected_genome_size', 'N/A')} bp\n")
                f.write(f"  Database name: {vc.get('database_name', 'N/A')}\n")
                f.write("\n")

            # Tool versions
            f.write("-" * 80 + "\n")
            f.write("TOOL VERSIONS\n")
            f.write("-" * 80 + "\n")
            for tool, version in sorted(data['tool_versions'].items()):
                f.write(f"  {tool}: {version}\n")
            f.write("\n")

            # Steps summary
            f.write("-" * 80 + "\n")
            f.write("PIPELINE STEPS\n")
            f.write("-" * 80 + "\n\n")

            completed = sum(1 for s in data['steps'] if s['status'] == 'completed')
            skipped = sum(1 for s in data['steps'] if s['status'] == 'skipped')
            failed = sum(1 for s in data['steps'] if s['status'] == 'failed')
            f.write(f"  Total steps: {len(data['steps'])}  "
                    f"(completed: {completed}, skipped: {skipped}, failed: {failed})\n\n")

            # List skipped steps prominently
            skipped_steps = [s for s in data['steps'] if s['status'] == 'skipped']
            if skipped_steps:
                f.write("  *** SKIPPED STEPS ***\n")
                for s in skipped_steps:
                    sample_str = f" [{s['sample']}]" if s.get('sample') else ""
                    f.write(f"    - {s['step']}{sample_str}: {s.get('notes', 'no reason given')}\n")
                f.write("\n")

            # Detailed step log
            for i, step in enumerate(data['steps'], 1):
                status_marker = {'completed': '[OK]', 'skipped': '[SKIP]', 'failed': '[FAIL]'}
                marker = status_marker.get(step['status'], f"[{step['status'].upper()}]")
                sample_str = f" [{step['sample']}]" if step.get('sample') else ""

                f.write(f"  Step {i}: {marker} {step['step']}{sample_str}\n")
                f.write(f"    Tool: {step['tool']}\n")
                f.write(f"    Time: {step['timestamp']}\n")
                if step.get('notes'):
                    f.write(f"    Note: {step['notes']}\n")
                if step['parameters']:
                    f.write(f"    Parameters:\n")
                    for pk, pv in step['parameters'].items():
                        f.write(f"      {pk}: {pv}\n")
                f.write("\n")

            f.write("=" * 80 + "\n")
            f.write("END OF PROVENANCE REPORT\n")
            f.write("=" * 80 + "\n")

        logging.info(f"Provenance report written to {path}")
        return path
