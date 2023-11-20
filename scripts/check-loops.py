def is_loop_exists(file_name):
    with open(file_name) as file:
        lines = [line.rstrip() for line in file]
        length = len(lines[1:])
        for i in range(1, length):
            line = lines[i]
            tokens = [int(token.strip()) for token in line.split()]
            tokens = tokens[1:]
            for token in tokens:
                if (i-1) == token:
                    print((i-1), tokens)
                    return True
    return False


def check_loops(root, file_name, name):
    n_loops = 0
    map = {}
    new_file_name = f"{root}indicators/{name}_data1.txt"
    new_file = open(new_file_name, "a")
    with open(file_name) as file:
        lines = [line.rstrip() for line in file]
        for line in lines:
            tokens = [token.strip() for token in line.split(",")]
            pattern = tokens[0]
            target = tokens[1]
            pattern_loop_exists = False
            target_loop_exists = False
            if pattern not in map:
                pattern_loop_exists = is_loop_exists(root+pattern)
                if pattern_loop_exists:
                    n_loops += 1
                map[pattern] = pattern_loop_exists
            else:
                pattern_loop_exists = map[pattern]
            if target not in map:
                target_loop_exists = is_loop_exists(root+target)
                if target_loop_exists:
                    n_loops += 1
                map[target] = target_loop_exists
            else:
                target_loop_exists = map[target]
            if (not pattern_loop_exists) and (not target_loop_exists):
                new_file.write(f"{pattern}, {target}\n")
    for m in map:
        if map[m] == True:
            print(m)

    return n_loops
    

name = "LV"
filename = (f"/Users/buddhi/Documents/Work/PhD/Implementation/mcsplit-plus-original/mcsplit-plus/indicators/{name}_data.txt")
root = "/Users/buddhi/Documents/Work/PhD/Implementation/mcsplit-plus-original/mcsplit-plus/"
print(check_loops(root, filename, name))
# print(is_loop_exists(root+"data/largerLV/targets/g58"))