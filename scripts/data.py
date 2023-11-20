def move(name, root):
    read_file_name = f"{root}/{name}_min_max_10.txt"
    write_file_name = f"{root}/{name}_data.txt"

    write_file = open(write_file_name, "a")
    with open(read_file_name) as file:
        lines = [line.rstrip() for line in file]
        for line in lines[1:]:
            tokens = [token.strip() for token in line.split(",")]
            write_file.write(f"{tokens[0]}, {tokens[1]}\n")


root = "/Users/buddhi/Documents/Work/PhD/Implementation/mcsplit-plus-original/mcsplit-plus/indicators"
files = ["largerLV", "si", "scalefree", "meshes-CVIU11", "LV", "images-CVIU11", "images-PR15"]
for file in files:
    move(file, root)
