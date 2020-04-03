import sys

output_bool = False

if len(sys.argv) == 1:
    print("Usage: combine_times <input>")
    print("Argument: input - the name of the input file containing run times.")
    print("Optional: output - an optional second argument is used to specify the output file. If output is not " +
          "specified, the program will only print results to the terminal.")
    print("Example usage without output: combine_times.py timing_file.txt")
    print("Example usage with output: combine_times.py timing_file.txt combined_timing.txt")

elif len(sys.argv) <= 3:
    input_file = sys.argv[1]
    if len(sys.argv) == 3:
        output_bool = True
    output_file = None
    if output_bool:
        output_file = sys.argv[3]

    real_seconds = 0.0
    user_seconds = 0.0
    sys_seconds = 0.0

    with open(input_file, "r") as input:
        for line in input:
            if "real" in line:
                time_list = line[line.find("real"):].split(" ")
                real = time_list[1][:-1].split("m")
                real_seconds += float(real[1]) + float(real[0]) * 60
                user = time_list[3][:-1].split("m")
                user_seconds += float(user[1]) + float(user[0]) * 60
                sys = time_list[5][:-2].split("m")
                sys_seconds += float(sys[1]) + float(sys[0]) * 60

    if output_bool:
        with open(output_file, "a") as output:
            output.write("real " + str(int(real_seconds)/60) + "m" + str(real_seconds % 60) + "s\n")
            output.write("user " + str(int(user_seconds) / 60) + "m" + str(user_seconds % 60) + "s\n")
            output.write("sys " + str(int(sys_seconds) / 60) + "m" + str(sys_seconds % 60) + "s\n")

    print("real " + str(int(real_seconds)/60) + "m" + str(real_seconds % 60) + "s\n" +
          "user " + str(int(user_seconds) / 60) + "m" + str(user_seconds % 60) + "s\n" +
          "sys " + str(int(sys_seconds) / 60) + "m" + str(sys_seconds % 60) + "s\n")
