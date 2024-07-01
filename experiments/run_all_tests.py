import subprocess

# Define the script to run and the different sets of arguments
script = 'comprehensive_test.py'
filenames = ['retention0s.csv', 'retention1s_1_dominated.csv', 'retention1s_2_dominated.csv', 'retention1s.csv', 'retention1s2.csv', 'retention1smerge.csv']

# Iterate over the arguments and run the script with each argument using -f
for filename in filenames:
    # Build the command
    command = ['python', script, '-f', filename]
    # Run the command
    result = subprocess.run(command, capture_output=True, text=True)
    # Print the output
    print(result.stdout)
    print(result.stderr)  # In case there are errors
