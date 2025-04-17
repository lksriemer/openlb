import os
import sys
import subprocess

# Path to your compiled C++ program
cmd = ["mpirun", "-np", "8", "--allow-run-as-root", "cylinder2d"]

def execute_simulation(order, nq):
    # Create the simulation output directory
    outdir = f"./simulation_order{order}_nq{nq}/"

    if ' ' in outdir:
        sys.exit(f"No spaces in '{outdir}' allowed")

    if not os.path.exists(outdir):
        os.mkdir(outdir)
    else:
        print(f"Directory {outdir} already exists. Output will be appended.")

    # Paths for logs
    stdout_log = os.path.join(outdir, "stdout.txt")
    stderr_log = os.path.join(outdir, "stderr.txt")

    # Construct the command
    simulation_cmd = cmd + [str(order), str(nq)]
    print("Executing:", ' '.join(simulation_cmd))

    # Open the log files
    with open(stdout_log, 'w') as stdout_file, open(stderr_log, 'w') as stderr_file:
        # Execute the command and write output to log files
        result = subprocess.run(simulation_cmd, stdout=stdout_file, stderr=stderr_file, text=True)

    # Read the last line of stdout log to extract drag coefficients
    extract_and_save_drag_coefficients(order, nq, stdout_log)

def extract_and_save_drag_coefficients(order, nq, stdout_log):
    # Read the stdout log file
    with open(stdout_log, 'r') as f:
        lines = f.readlines()

    if not lines:
        print(f"No output found in {stdout_log} for order={order}, nq={nq}")
        return

    # Get the last non-empty line
    for line in reversed(lines):
        if line.strip():
            last_line = line.strip()
            break
    else:
        print(f"No non-empty lines found in {stdout_log} for order={order}, nq={nq}")
        return

    if 'Drag coefficients:' in last_line:
        # Parse the last line
        parts = last_line.strip().split(',')
        if len(parts) >= 2:
            mean_part = parts[0]
            std_part = parts[1]

            mean_drag = mean_part.split(':')[-1].strip()
            std_drag = std_part.split(':')[-1].strip()

            # Save the results to a summary file
            # summary_file = "drag_coefficients.txt"
            file_exists = os.path.isfile(summary_file)
            with open(summary_file, 'a') as f_out:
                if not file_exists:
                    f_out.write('Order\tNq\tMean_Drag\tStd_Drag\n')
                f_out.write(f"{order}\t{nq}\t{mean_drag}\t{std_drag}\n")
            print(f"Saved drag coefficients for order={order}, nq={nq}")
        else:
            print(f"Unexpected format in last line: {last_line}")
    else:
        print(f"Drag coefficients not found in last line for order={order}, nq={nq}")

if __name__ == '__main__':

    # Remove the existing summary file if you want to start fresh
    summary_file = "drag_coefficients.txt"
    if os.path.exists(summary_file):
        os.remove(summary_file)

    # simulations = [
    #     (1, 2),
    #     (2, 10),
    #     (3, 20),
    #     (4, 40),
    #     (5, 80)
    # ]

        simulations = [
        (1, 3),
        (2, 5),
        (3, 7),
        (4, 9),
        (5, 11),
        (6, 13),
        (7, 15),
        (8, 17),
        (9, 19),
        (10, 21),
        (11, 23),
        (12, 25),
        (13, 27),
        (14, 29)
    ]

    for order, nq in simulations:
        execute_simulation(order, nq)