import argparse
import os
import requests
import subprocess

def create_file(out_filename,in_url):
    r = requests.get(in_url, allow_redirects=True)
    open(out_filename, 'wb').write(r.content)

def run_metric(output_dir, name, case_pos, ctrl_pos):
    # Create the output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, f'{name}.log.txt')

    # Download R script to benchmark folder
    R_script_url = "https://raw.githubusercontent.com/sorensandgaard/ob_anonymization_R3_clustering/main/clustering.R"
    script_R_file = os.path.join(output_dir, f'clustering.R')
    create_file(script_R_file,R_script_url)

    # Run R script
    wrapper_R = f"envs/R_wrapper.sh"
    outfile_pos = f"{output_dir}/{name}.somefile.json"
    R_command = f"{wrapper_R} {script_R_file} {case_pos} {ctrl_pos} {outfile_pos}"
    content = f"R Command\n{R_command}\n\n"
    a = subprocess.run(R_command.split(),capture_output=True,text=True)
    content += f"RCD output:\n{a.stdout}\n"

    # cleanup_command = f"rm {script_R_file}"
    # a = subprocess.run(cleanup_command.split(),capture_output=True,text=True)

    with open(log_file, 'w') as file:
        file.write(content)


def main():
    # Create argument parser
    parser = argparse.ArgumentParser(description='Run metic on files.')

    # Add arguments
    parser.add_argument('--output_dir', type=str, help='output directory where metic will store results.')
    parser.add_argument('--name', type=str, help='name of the dataset')
    parser.add_argument('--case_expr', type=str, help='Anonymized version of data')
    parser.add_argument('--ctrl_expr', type=str, help='Control version of data')

    # Parse arguments
    args = parser.parse_args()

    case_pos = getattr(args, 'case_expr')
    ctrl_pos = getattr(args, 'ctrl_expr')

    run_metric(args.output_dir, args.name, case_pos, ctrl_pos)


if __name__ == "__main__":
    main()
