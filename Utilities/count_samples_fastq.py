import sys

output_bool = False

if len(sys.argv) == 1:
    print("Usage: count_sample_fastq <input>")
    print("Argument: input - the name of the input file in fastq format.")
    print("Optional: output - an optional second argument is used to specify the output file. If output is not " +
          "specified, the program will only print results to the terminal.")
    print("Example usage without output: count_sample_fastq.py input.fastq")
    print("Example usage with output: count_sample_fastq.py input.fastq output.txt")

elif len(sys.argv) <= 3:
    input_file = sys.argv[1]
    if len(sys.argv) == 3:
        output_bool = True
    output_file = None
    if output_bool:
        output_file = sys.argv[2]

    counter = 0

    with open(input_file, "r") as input:
        for line in input:
            counter += 1

    if output_bool:
        with open(output_file, "a") as output:
            output.write(str(counter/4) + " samples in " + input_file + ".")

    print(str(counter/4) + " samples in " + input_file + ".")
