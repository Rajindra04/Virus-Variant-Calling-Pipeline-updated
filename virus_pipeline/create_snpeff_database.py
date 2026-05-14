import sys
import shutil
import argparse
import os
import subprocess
import logging

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

def run_command(command):
    logging.info(f"Running command: {command}")
    process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    stdout, stderr = process.communicate()
    if process.returncode != 0:
        raise Exception(f"Command execution failed with return code {process.returncode}, stderr: {stderr.decode('utf-8')}")
    return stdout.decode("utf-8"), stderr.decode("utf-8")

def validate_genbank_file(genbank_file):
    try:
        with open(genbank_file, "r") as f:
            content = f.read()
            if not content.startswith("LOCUS"):
                raise ValueError("GenBank file does not start with 'LOCUS'. Invalid format.")
            if "FEATURES" not in content:
                raise ValueError("GenBank file lacks 'FEATURES' section. Annotations may be missing.")
        logging.info(f"GenBank file {genbank_file} appears valid.")
    except Exception as e:
        raise Exception(f"GenBank file validation failed: {str(e)}")

def validate_fasta_genbank_match(genbank_file, reference_fasta):
    try:
        with open(reference_fasta, "r") as f:
            fasta_id = f.readline().strip().lstrip(">").split()[0]
        
        with open(genbank_file, "r") as f:
            for line in f:
                if line.startswith("LOCUS"):
                    gbk_id = line.split()[1]
                    break
            else:
                raise ValueError("No LOCUS line found in GenBank file.")
        
        if fasta_id != gbk_id:
            raise ValueError(f"Sequence ID mismatch: FASTA ID '{fasta_id}' does not match GenBank ID '{gbk_id}'.")
        logging.info(f"Sequence IDs match: {fasta_id}")
    except Exception as e:
        raise Exception(f"FASTA-GenBank validation failed: {str(e)}")

def create_snpeff_config(output_dir, database_name, reference_fasta):
    data_dir = os.path.abspath(os.path.join(output_dir, "data"))
    custom_config = os.path.join(output_dir, "snpEff.config")
    config_content = f"""# SnpEff configuration file
data.dir = {data_dir}

# Database entry
{database_name}.genome = {database_name}
{database_name}.reference = {os.path.abspath(reference_fasta)}
"""
    with open(custom_config, "w") as config_file:
        config_file.write(config_content)
    
    logging.info(f"SnpEff configuration file created: {custom_config}")
    return custom_config

def prepare_files(genbank_file, reference_fasta, output_dir, database_name):
    data_dir = os.path.join(output_dir, "data", database_name)
    os.makedirs(data_dir, exist_ok=True)
    fasta_dest = os.path.join(data_dir, "sequences.fa")
    genbank_dest = os.path.join(data_dir, "genes.gbk")
    run_command(f"cp {reference_fasta} {fasta_dest}")
    run_command(f"cp {genbank_file} {genbank_dest}")
    if not os.path.exists(fasta_dest):
        raise FileNotFoundError(f"FASTA file not found at {fasta_dest}")
    if not os.path.exists(genbank_dest):
        raise FileNotFoundError(f"GenBank file not found at {genbank_dest}")
    logging.info(f"Files prepared: {fasta_dest}, {genbank_dest}")

def build_snpeff_database(database_name, config_file, output_dir, java_path):
    data_dir = os.path.abspath(os.path.join(output_dir, "data"))
    
    # 1. Dynamically locate the JAR using the active Python environment prefix
    # In Conda, this is usually: envs/dengue_pipeline/share/snpeff/snpEff.jar
    conda_prefix = sys.prefix
    potential_jar = os.path.join(conda_prefix, "share", "snpeff", "snpEff.jar")
    
    # Fallback search if the versioned folder is used (e.g., snpeff-5.2-0)
    if not os.path.exists(potential_jar):
        share_dir = os.path.join(conda_prefix, "share")
        if os.path.exists(share_dir):
            for folder in os.listdir(share_dir):
                if folder.startswith("snpeff"):
                    test_path = os.path.join(share_dir, folder, "snpEff.jar")
                    if os.path.exists(test_path):
                        potential_jar = test_path
                        break

    # 2. Verify we actually found it
    if not os.path.exists(potential_jar):
        logging.error(f"Could not find snpEff.jar at {potential_jar}")
        # Last ditch effort: check if it's in the current directory
        if os.path.exists("snpEff.jar"):
            potential_jar = "snpEff.jar"
        else:
            raise FileNotFoundError("snpEff.jar not found in Conda share/ or current directory.")

    # 3. Construct the command using the absolute path to the JAR
    # Use quotes around the jar path in case there are spaces
    build_command = (
        f"{java_path} -Xmx4g -jar '{potential_jar}' build -genbank -v {database_name} "
        f"-noCheckCds -noCheckProtein "
        f"-c {config_file} "
        f"-dataDir {data_dir}"
    )
    
    stdout, stderr = run_command(build_command)
    logging.info("SnpEff build output: %s", stdout)
    logging.info("SnpEff build error: %s", stderr)

def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]
    
    parser = argparse.ArgumentParser()
    parser.add_argument("--genbank_file", type=str, required=True, help="Path to GenBank (.gb or .gbk) file.")
    parser.add_argument("--reference_fasta", type=str, required=True, help="Path to reference FASTA file.")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to output directory for SnpEff database.")
    parser.add_argument("--database_name", type=str, default="denv1", help="Name of the SnpEff database.")
    parser.add_argument("--skip_id_validation", action="store_true", help="Skip validation of sequence ID matching.")
    
    # FIX: Added the missing --snpeff_java argument
    parser.add_argument("--snpeff_java", type=str, default="java", help="Path to specific Java version for SnpEff.")
    
    args = parser.parse_args(argv)

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    try:
        validate_genbank_file(args.genbank_file)
        if not args.skip_id_validation:
            validate_fasta_genbank_match(args.genbank_file, args.reference_fasta)
            
        config_file = create_snpeff_config(args.output_dir, args.database_name, args.reference_fasta)
        prepare_files(args.genbank_file, args.reference_fasta, args.output_dir, args.database_name)
        
        # FIX: Pass the java path to the build function
        build_snpeff_database(args.database_name, config_file, args.output_dir, args.snpeff_java)
        
        logging.info(f"SnpEff database '{args.database_name}' created successfully in {args.output_dir}")
    except Exception as e:
        logging.error(f"Error occurred during database creation: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()
